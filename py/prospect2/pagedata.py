#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.pagedata
==================

PageData: the data of one prospect2 page, a lightweight wrapper in the style of
desispec Spectra, with the same filter / stack / write workflow::

    from prospect2 import read_page, write_page, stack
    pages = [read_page(f) for f in files]
    pages = [p[p.redshifts['SPECTYPE'] == 'QSO'] for p in pages]
    write_page('qso.html', stack(pages))

Pages made from FITS go through the same object: make_html(spectra, ...) is
PageData.from_spectra(spectra, ...).to_html().

Per-pixel data stay in their quantized integer form ("q2", see
prospect2.encode), so selecting and stacking are lossless: a page stacked from
pieces has the same arrays as one generated directly from the same spectra.
Per-target metadata are astropy Tables (fibermap, redshifts) with DESI column
names, plus `extra` for prospect2's own per-target values.  Display values
(target-bit names, normalized exposure time, ...) are not stored here; they
are recomputed from the tables whenever a page is written.
"""

import re, json, hashlib, warnings
import numpy as np
from astropy.table import Table, Column, MaskedColumn, vstack

from desiutil.log import get_logger

from ._version import __version__
from .encode import quantize, quantize_model, pack, unpack, ARRAY_INFO, decode_noise_flux, decode_model
from .models import camera_models, merge_camera_models, alt_fit_models
from .grouping import grouping_columns, concat_columns, key_columns
from .util import fibermap_value
from .io import read_lines

#- redrock's "same redshift" criterion for marking near-duplicate alternate fits
_DUP_KMS = 1000

#- stored per-target columns (only the ones the viewer uses, plus raw target bits)
REDSHIFT_COLUMNS = ('Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'DELTACHI2', 'CHI2', 'NPIXELS')
#- grouping values: from fibermap columns or headers (grouping.grouping_columns)
GROUP_COLUMNS = {'SURVEY': str, 'PROGRAM': str, 'TILEID': np.int32, 'FIBER': np.int32,
                 'NIGHT': np.int32, 'EXPID': np.int32}
FIBERMAP_COLUMNS = ('TARGET_RA', 'TARGET_DEC', 'OBJTYPE', 'MORPHTYPE',
                    'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_R', 'EBV',
                    'COADD_EXPTIME', 'EXPTIME', 'COADD_NUMEXP', 'COADD_FIBERSTATUS', 'FIBERSTATUS')
TARGET_BIT_COLUMN = re.compile(r'^((SV[123]_)?(DESI|BGS|MWS|SCND)|CMX)_TARGET$')

#- lowercase per-target keys: stored (in PageData.extra) or derived (recomputed when written)
STORED_KEYS = ('snr', 'masked_cams', 'zfits', 'src')
DERIVED_KEYS = ('classes', 'primary', 'zwarn_names', 'fiberstatus', 'fiberstatus_names', 'exptime', 'numexp')
SOURCE_KEYS = ('source', 'redrock', 'model_source', 'rrdetails')

FLUX_UNIT = '1e-17 erg/s/cm²/Å'


# ---------------------------------------------------------------- payload format version

#- The payload format is versioned MAJOR.MINOR (doc/implementation.md "Payload versioning"):
#- major = a stored key removed/renamed or changed in meaning or encoding; minor = stored keys
#- added.  PageData reads any page of the same major; from a newer minor it drops, with a warning,
#- the keys it doesn't know (user decision 2026-11-18, "middle option").  The viewer template
#- must display every page of its major (user decision 2026-11-18): it tolerates missing
#- minor-added keys, and anything that needs them requires a major bump.
#- test_pagedata.py::test_payload_keyset_guard enforces bumps.
FORMAT_VERSION = '1.0'


def parse_version(v):
    """'MAJOR.MINOR' -> (major, minor) ints; ValueError if malformed"""
    try:
        major, minor = str(v).split('.')
        return int(major), int(minor)
    except (ValueError, TypeError):
        raise ValueError(f'malformed format_version {v!r} (expected "MAJOR.MINOR")') from None


def check_readable(v):
    """
    Raise ValueError unless this code can read (load, select, stack) a payload
    of format version `v`: same major as FORMAT_VERSION.  Returns True if `v`
    is a newer minor, whose added keys this code doesn't know (from_payload
    drops them with a warning), else False.
    """
    if v is None:
        raise ValueError('page has no format_version (made by a prototype before format 1.0); regenerate it')
    major, minor = parse_version(v)
    kmajor, kminor = parse_version(FORMAT_VERSION)
    if major != kmajor:
        raise ValueError(f'page data format {v} is not readable by this prospect2 (reads {kmajor}.x); '
                         + ('upgrade prospect2' if major > kmajor else 'regenerate the page'))
    return minor > kminor


#- keys of the current format that from_payload reads; anything else in a newer-minor page is dropped
TOP_KEYS = frozenset(('format_version', 'generator', 'title', 'vi_key', 'sources', 'spgrp', 'keycols', 'wave0',
                      'dwave', 'nwave', 'nspec', 'flux_unit', 'columns', 'targets', 'lines', 'overlaps',
                      'alt_models', 'arrays'))
TABLES = ('fibermap', 'redshifts')
ARRAY_KEYS = frozenset(('rows', 'pack', 'enc', 'data', 'scale', 'owner'))


def _drop_newer(d):
    """
    Copy of payload `d` (a newer minor of this major) without the keys this
    code doesn't know, and a list of what was dropped (unknown lowercase
    per-target keys may be stored or derived; this code can't tell, and
    ignores both).
    """
    dropped = [k for k in d if k not in TOP_KEYS]
    d = {k: v for k, v in d.items() if k in TOP_KEYS}
    columns = dict(d['columns'])
    targets = dict(d['targets'])
    for t in [t for t in columns if t not in TABLES]:
        dropped.append(f'columns table {t}')
        for c in columns.pop(t) or {}:
            targets.pop(c, None)
    d['columns'], d['targets'] = columns, targets
    dropped += [f'per-target {k}' for k in targets
                if k != k.upper() and k not in STORED_KEYS and k not in DERIVED_KEYS]
    extra = sorted({k for s in d['sources'] for k in s if k not in SOURCE_KEYS})
    dropped += [f'sources field {k}' for k in extra]
    d['sources'] = [{k: v for k, v in s.items() if k in SOURCE_KEYS} for s in d['sources']]
    arrays = {}
    for name, a in d['arrays'].items():
        if name not in ARRAY_INFO:
            dropped.append(f'array {name}')
            continue
        dropped += [f'array {name} field {k}' for k in a if k not in ARRAY_KEYS]
        arrays[name] = {k: v for k, v in a.items() if k in ARRAY_KEYS}
    d['arrays'] = arrays
    return d, dropped


# ---------------------------------------------------------------- helpers

def _num(x, ndigits):
    """round(float(x), ndigits), or None if x is None or not finite"""
    if x is None or not np.isfinite(float(x)):
        return None
    return round(float(x), ndigits)


def target_classes(fibermap, i):
    """(all target-class names, primary classes for badges) of fibermap row i, from the main-survey bits"""
    try:
        from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
    except ImportError:
        return [], []
    names = []
    for col, mask in [('DESI_TARGET', desi_mask), ('BGS_TARGET', bgs_mask),
                      ('MWS_TARGET', mws_mask), ('SCND_TARGET', scnd_mask)]:
        bits = int(fibermap_value(fibermap, i, col, default=0))
        if bits != 0:
            names += [n for n in mask.names(bits)
                      if n not in ('BGS_ANY', 'MWS_ANY', 'SCND_ANY') and not n.endswith(('_SOUTH', '_NORTH'))]
    #- primary classes for compact badges
    primary = [n for n in names if n in ('ELG', 'LRG', 'QSO', 'BGS_BRIGHT', 'BGS_FAINT', 'MWS_BROAD',
                                         'MWS_MAIN_BLUE', 'MWS_MAIN_RED', 'STD_FAINT', 'STD_BRIGHT',
                                         'STD_WD', 'MWS_WD', 'ELG_LOP', 'ELG_VLO', 'ELG_HIP')]
    return names, primary


def _fit_table(d):
    """
    Display list for one target's raw rrdetails records `d`: [[z, spectype,
    subtype, chi2 - best chi2, zwarn, dup], ...] (compact arrays to keep the page
    small; expanded to objects in JS).  dup = index of a better-ranked fit within
    1000 km/s, else -1.
    """
    C = 299792.458
    s = lambda x: (x.decode() if isinstance(x, bytes) else str(x)).strip()
    out = []
    for k in range(len(d)):
        zk = float(d['z'][k])
        dup = next((j for j in range(k) if C * abs(zk - d['z'][j]) / (1 + d['z'][j]) < _DUP_KMS), -1)
        out.append([round(zk, 6), s(d['spectype'][k]), s(d['subtype'][k]).strip('-'),
                    round(float(d['chi2'][k] - d['chi2'][0]), 1), int(d['zwarn'][k]), dup])
    return out


def _template_versions(zcat):
    """{TEMNAMnn/TEMVERnn: value} from a redrock table's meta"""
    return {k: v for k, v in zcat.meta.items() if re.fullmatch(r'TEM(NAM|VER)\d+', k)}


def _combine(spectra, zcat, models, zfits):
    """
    Stack lists of inputs (one entry per input file) into single objects,
    raising ValueError if the redrock template versions differ (models computed
    from templates must use the templates of the fit).  Also returns the
    per-row grouping values (grouping.grouping_columns), read before stacking
    because stack keeps only the first input's headers.
    """
    from desispec.spectra import stack as stack_spectra
    n = len(spectra)
    for name, x in (('zcat', zcat), ('models', models), ('zfits', zfits)):
        if x is not None and (not isinstance(x, (list, tuple)) or len(x) != n):
            raise ValueError(f'spectra is a list of {n}, so {name} must be a list of {n} too')
    zcats = zcat if zcat is not None else [sp.redshifts for sp in spectra]
    if all(z is not None for z in zcats):
        versions = [_template_versions(z) for z in zcats]
        for i, v in enumerate(versions[1:], 1):
            if v != versions[0]:
                diff = {k: (versions[0].get(k), v.get(k)) for k in sorted(set(v) | set(versions[0]))
                        if versions[0].get(k) != v.get(k)}
                raise ValueError(f'redrock template versions differ between inputs 0 and {i} '
                                 f'(keyword: (input 0, input {i})): {diff}')
    elif any(z is not None for z in zcats):
        raise ValueError('some inputs have redshifts and some do not')
    groups = concat_columns([grouping_columns(s) for s in spectra])
    sp = stack_spectra(list(spectra))
    if zcat is not None:
        zcat = vstack(list(zcat), metadata_conflicts='silent')
    if models is not None:
        models = {b: np.concatenate([m[b] for m in models]) for b in models[0]}
    if zfits is not None:
        zfits = [f for flist in zfits for f in flist]
    return sp, zcat, models, zfits, groups


def _objcol(name, values):
    """astropy Column of Python objects (lists, dicts, None), one per row"""
    a = np.empty(len(values), dtype=object)
    for i, v in enumerate(values):
        a[i] = v
    return Column(a, name=name)


def _clean(col):
    """Copy of a 1-d table column, bytes decoded and strings stripped, masked only if anything is"""
    if col.ndim != 1:
        raise ValueError(f'column {col.name} is not 1-d')
    data = np.ma.getdata(col)
    mask = np.ma.getmaskarray(col)
    if data.dtype.kind == 'S':
        data = np.char.decode(data, 'utf-8')
    if data.dtype.kind == 'U':
        data = np.char.strip(data)
    else:
        data = np.array(data)
    if mask.any():
        return MaskedColumn(data, name=col.name, mask=mask)
    return Column(data, name=col.name)


def _column_from_values(name, values, dtype):
    """Column from Python values with None = masked; None if every value is None"""
    mask = np.array([v is None for v in values], dtype=bool)
    if mask.all():
        return None
    fill = '' if dtype is str else 0
    data = np.array([fill if v is None else v for v in values], dtype=dtype)
    return MaskedColumn(data, name=name, mask=mask) if mask.any() else Column(data, name=name)


def _always_string(name):
    """int64 columns written as JSON strings regardless of value (int64 > 2**53 is not exact in JS)"""
    return name == 'TARGETID' or TARGET_BIT_COLUMN.match(name) is not None


def _json_column(col):
    """(JSON list, dtype code) for a 1-d table column; masked and non-finite values -> null"""
    name = col.name
    if col.ndim != 1:
        raise ValueError(f'column {name} is not 1-d; only 1-d columns can be written to a page')
    data = np.ma.getdata(col)
    mask = np.ma.getmaskarray(col)
    kind = data.dtype.kind
    if kind in 'SU':
        if kind == 'S':
            data = np.char.decode(data, 'utf-8')
        vals, code = [str(v) for v in data], f'U{max(1, data.dtype.itemsize // 4)}'   #- keeps the width
    elif kind == 'b':
        vals, code = [bool(v) for v in data], 'bool'
    elif kind in 'iu':
        code = f'{kind}{data.dtype.itemsize}'
        big = data.dtype.itemsize == 8 and (_always_string(name) or bool(np.any(data > 2**53)) or bool(np.any(data < -2**53)))
        vals = [str(int(v)) for v in data] if big else [int(v) for v in data]
    elif kind == 'f':
        code = f'f{data.dtype.itemsize}'
        #- float32: shortest decimal that round-trips in float32 (not the float64 expansion)
        conv = (lambda v: float(str(np.float32(v)))) if data.dtype.itemsize == 4 else float
        vals = [conv(v) if np.isfinite(v) else None for v in data]
    else:
        raise ValueError(f'column {name}: dtype {data.dtype} can not be written to a page')
    return [None if m else v for v, m in zip(vals, mask)], code


def _column_from_json(name, values, code):
    """Inverse of _json_column: floats null -> NaN; other nulls -> masked"""
    if code.startswith('f'):
        return Column(np.array([np.nan if v is None else float(v) for v in values], dtype=code), name=name)
    if code.startswith('U'):
        dtype = np.dtype(code)
    elif code == 'bool':
        dtype = bool
    else:
        dtype = np.dtype(code)
        values = [None if v is None else int(v) for v in values]
    col = _column_from_values(name, values, dtype)
    if col is None:                         #- all null
        return MaskedColumn(np.zeros(len(values), dtype=dtype), name=name, mask=True)
    return col


def _take_rows(arr, rows):
    """Rows `rows` of a quantized array dict (i8, esc, per-row side arrays); escapes gathered in order"""
    i8 = arr['i8']
    offs = np.concatenate([[0], np.cumsum((i8 == -128).sum(axis=1))])
    esc = arr['esc']
    out = dict(arr)
    out['i8'] = i8[rows]
    out['esc'] = (np.concatenate([esc[offs[r]:offs[r + 1]] for r in rows]) if len(rows) else esc[:0]).astype('<i4')
    for k in _side_keys(arr):
        out[k] = np.asarray(arr[k])[rows]
    return out


def _side_keys(arr):
    """per-row side arrays of a quantized array dict (e.g. scale, owner)"""
    return [k for k in arr if k not in ('rows', 'enc', 'i8', 'esc')]


# ---------------------------------------------------------------- PageData

class PageData:
    """
    The data of one prospect2 page.  Make one with PageData.from_spectra(),
    PageData.from_html() / prospect2.io.read_page(), or PageData.from_payload().

    Attributes:
        fibermap: astropy Table, one row per spectrum, DESI column names (a subset
            of the input fibermap plus the grouping values SURVEY, PROGRAM,
            TILEID, FIBER, NIGHT, EXPID where known).  Editable; edits are written.
        redshifts: astropy Table (TARGETID + redrock columns) or None.  Editable.
        extra: astropy Table of prospect2 per-target values that can't be
            recomputed: snr ({band: median S/N}), masked_cams, zfits (display
            list of the top redrock fits, or None), src (index into sources).
        sources: list of dicts (source, redrock, model_source, rrdetails): the
            input files of the page(s) these spectra come from.
        arrays: quantized per-pixel arrays {name: dict(rows=, enc=, i8=, esc=, ...)}.
        title: the page title (the label at the top left of the page), or None
            for the default (the single source's file name, else 'prospect2';
            write_page uses the file name for several sources).  Settable;
            kept by selection, not by stack (see its title argument).
        wave0, dwave, nwave: the linear wavelength grid; spgrp, overlaps,
        alt_models, flux_unit.

    Decoded (lazily, float64, NaN where masked): wave, flux, noise, model.
    Index with a bool mask, integer indices or a slice: p[p.fibermap['FIBER'] < 2010].
    """

    def __init__(self, fibermap, redshifts, extra, arrays, wave0, dwave, nwave, sources,
                 spgrp=None, overlaps=(), alt_models=None, flux_unit=FLUX_UNIT, title=None):
        n = len(fibermap)
        if 'TARGETID' not in fibermap.colnames:
            raise ValueError('fibermap must have a TARGETID column')
        for name, tab in (('redshifts', redshifts), ('extra', extra)):
            if tab is not None and len(tab) != n:
                raise ValueError(f'{name} has {len(tab)} rows for {n} spectra')
        for name, a in arrays.items():
            nrow = len(a['i8'])
            if a['rows'] == 'spec' and nrow != n:
                raise ValueError(f'array {name} has {nrow} rows for {n} spectra')
            if a['rows'] != 'spec' and 'owner' not in a:
                raise ValueError(f'array {name} (rows {a["rows"]!r}) has no owner')
            if a['i8'].shape[1:] != (nwave,):
                raise ValueError(f'array {name} is not {nwave} pixels wide')
            if 'owner' in a:
                owner = _owner(a['owner'])
                if len(owner) and not (0 <= owner[:, 0].min() and owner[:, 0].max() < n):
                    raise ValueError(f'array {name}: owner spectrum index outside 0..{n - 1}')
                if a['rows'] == 'alt':        #- (spectrum, index into that spectrum's zfits)
                    zfits = extra['zfits'] if extra is not None and 'zfits' in extra.colnames else None
                    for i, k in owner[:, :2]:
                        if not (zfits is not None and zfits[i] is not None and 0 <= k < len(zfits[i])):
                            raise ValueError(f'array {name}: owner ({i}, {k}) is not a fit in zfits')
        if extra is not None and 'src' in extra.colnames and len(extra) and \
                not (0 <= np.min(extra['src']) and np.max(extra['src']) < len(sources)):
            raise ValueError(f'extra["src"] is not an index into the {len(sources)} sources')
        self.fibermap = fibermap
        self.redshifts = redshifts
        self.extra = extra
        self.arrays = arrays
        self.wave0, self.dwave, self.nwave = float(wave0), float(dwave), int(nwave)
        self.sources = [dict(s) for s in sources]      #- copies: selections don't share them
        self.spgrp = spgrp
        self.overlaps = [list(o) for o in overlaps]
        self.alt_models = dict(alt_models) if alt_models is not None else None
        self.flux_unit = flux_unit
        self.title = title
        self._decoded = {}
        #- the row order the arrays and extra are in, to catch tables reordered in place
        self._targetids = np.array(fibermap['TARGETID'], dtype=np.int64)

    def __len__(self):
        return len(self.fibermap)

    def num_spectra(self):
        return len(self)

    def __repr__(self):
        return (f'<PageData: {len(self)} spectra, {len(self.sources)} source(s), '
                f'arrays {", ".join(self.arrays)}>')

    # ------------------------------------------------------------ decoded arrays

    @property
    def wave(self):
        return self.wave0 + self.dwave * np.arange(self.nwave)

    def _decode(self):
        if 'flux' not in self._decoded:
            self._decoded['flux'], self._decoded['noise'] = decode_noise_flux(self.arrays['noise'], self.arrays['flux'])
        return self._decoded

    @property
    def flux(self):
        """[nspec, nwave] float64, NaN where masked"""
        return self._decode()['flux']

    @property
    def noise(self):
        """[nspec, nwave] float64 (1 sigma), NaN where masked"""
        return self._decode()['noise']

    @property
    def model(self):
        """[nspec, nwave] best-fit redrock model, or None"""
        if 'model' not in self.arrays:
            return None
        if 'model' not in self._decoded:
            self._decoded['model'] = decode_model(self.arrays['model'])
        return self._decoded['model']

    # ------------------------------------------------------------ selection

    def __getitem__(self, sel):
        self._check_rows()      #- else a table sorted in place would become the new row order
        n = len(self)
        if isinstance(sel, np.ma.MaskedArray):
            sel = sel.filled(False)
        if not isinstance(sel, slice):
            sel = np.asarray(sel)
            if sel.size == 0:
                sel = sel.astype(int)              #- p[[]]: an empty list is float to numpy
        idx = np.arange(n)[sel]
        idx = np.atleast_1d(idx)
        arrays = {}
        for name, a in self.arrays.items():
            if a['rows'] == 'spec':
                arrays[name] = _take_rows(a, idx)
            else:
                #- rows owned by spectra (owner column 0): keep those of kept spectra, in new order
                owner = _owner(a['owner'])
                byspec = {}
                for r, i in enumerate(owner[:, 0]):
                    byspec.setdefault(int(i), []).append(r)
                rows = [r for i in idx for r in byspec.get(int(i), [])]
                newi = [j for j, i in enumerate(idx) for _ in byspec.get(int(i), [])]
                arrays[name] = _take_rows(a, np.array(rows, dtype=int))
                newowner = owner[np.array(rows, dtype=int)].copy()
                newowner[:, 0] = newi
                arrays[name]['owner'] = newowner
        extra = self.extra[idx]
        #- keep only the sources still used, in their original order
        used = sorted(set(int(s) for s in extra['src'])) if 'src' in extra.colnames else [0]
        remap = {s: j for j, s in enumerate(used)}
        if 'src' in extra.colnames:
            extra['src'] = [remap[int(s)] for s in extra['src']]
        return PageData(self.fibermap[idx], None if self.redshifts is None else self.redshifts[idx], extra, arrays,
                        self.wave0, self.dwave, self.nwave, [self.sources[s] for s in used],
                        spgrp=self.spgrp, overlaps=self.overlaps, alt_models=self.alt_models,
                        flux_unit=self.flux_unit, title=self.title)

    # ------------------------------------------------------------ construction

    @classmethod
    def from_spectra(cls, spectra, zcat=None, models=None, zfits=None, title=None, provenance=None,
                     alt_models=2, alt_max_dchi2=100, model_precision=100, alt_model_precision=30,
                     groups=None):
        """
        PageData from desispec Spectra; see prospect2.page.make_html for the
        arguments.  `spectra` may be a list (one entry per input file), with
        zcat/models/zfits lists of the same length.  `groups` (internal) are the
        per-row grouping values of stacked inputs; default
        grouping.grouping_columns(spectra).
        """
        log = get_logger()
        if isinstance(spectra, (list, tuple)):
            spectra, zcat, models, zfits, groups = _combine(spectra, zcat, models, zfits)
        sp = spectra
        provenance = provenance or {}
        if sp.fibermap is None or 'TARGETID' not in sp.fibermap.colnames:
            raise ValueError('spectra.fibermap must have a TARGETID column (the only required column)')
        nspec = sp.num_spectra()
        if nspec == 0:
            raise ValueError('no spectra')

        #- explicit arguments win over what's attached to spectra
        if zcat is None:
            zcat = sp.redshifts
        elif sp.redshifts is not None:
            log.info('zcat argument overrides spectra.redshifts')
        if models is None:
            models = sp.model
        elif sp.model is not None:
            log.info('models argument overrides spectra.model')

        #- inputs must be row-matched to spectra
        tids = np.asarray(sp.fibermap['TARGETID'])
        if groups is None:
            groups = grouping_columns(sp)
        spgrps = sorted({g for g in groups['SPGRP'] if g is not None})
        if len(spgrps) > 1:
            raise ValueError(f'inputs mix spectral groupings (SPGRP): {", ".join(spgrps)}')
        spgrp = spgrps[0] if spgrps else None
        if zcat is not None:
            ztids = np.asarray(zcat['TARGETID'])
            if len(ztids) != nspec or np.any(ztids != tids):
                bad = next((i for i in range(min(len(ztids), nspec)) if ztids[i] != tids[i]), min(len(ztids), nspec))
                raise ValueError(f'zcat is not row-matched to spectra: {len(ztids)} vs {nspec} rows, first mismatch at row {bad}')
        if models is not None:
            if zcat is None:
                raise ValueError('models given without zcat')
            for b in sp.bands:
                if b not in models or np.shape(models[b]) != np.shape(sp.flux[b]):
                    raise ValueError(f'models[{b!r}] does not match spectra.flux[{b!r}] in shape')
        if zfits is not None:
            if zcat is None:
                raise ValueError('zfits given without zcat')
            if len(zfits) != nspec:
                raise ValueError(f'zfits has {len(zfits)} entries for {nspec} spectra')

        #- Per-camera median S/N (before merging cameras)
        snr = {}
        for b in sp.bands:
            with np.errstate(invalid='ignore'):
                s = sp.flux[b] * np.sqrt(sp.ivar[b])
            with warnings.catch_warnings():         #- fully masked fibers -> NaN (shown as "—")
                warnings.simplefilter('ignore', RuntimeWarning)
                snr[b] = np.nanmedian(np.where(sp.ivar[b] > 0, s, np.nan), axis=1)

        #- Camera overlap regions [lo, hi, label], from the per-camera grids
        bands = sorted(sp.bands, key=lambda b: sp.wave[b][0])
        overlaps = [[float(sp.wave[b2][0]), float(sp.wave[b1][-1]), f'{b1}/{b2} overlap']
                    for b1, b2 in zip(bands[:-1], bands[1:]) if sp.wave[b2][0] < sp.wave[b1][-1]]

        #- Redrock model on per-camera grids, then merged
        model = None
        model_source = None
        if zcat is not None:
            if models is not None:
                cam_models = models
                model_source = provenance.get('model_source', 'input models')
            else:
                cam_models = camera_models(sp, zcat)
                model_source = 'redrock templates + COEFF'
                log.info('redrock models computed from templates + COEFF')

        #- Merge cameras onto a single linear grid
        from desispec.coaddition import coadd_cameras
        sp_all = coadd_cameras(sp)
        band = sp_all.bands[0]
        wave = sp_all.wave[band]
        flux = sp_all.flux[band].copy()
        ivar = sp_all.ivar[band]
        bad = (ivar <= 0)
        if sp_all.mask is not None:
            bad |= (sp_all.mask[band] != 0)
        flux[bad] = np.nan
        with np.errstate(divide='ignore'):
            noise = np.where(bad, np.nan, 1 / np.sqrt(ivar))
        if zcat is not None:
            model = merge_camera_models(sp, cam_models, wave)

        dw = np.diff(wave)
        assert np.allclose(dw, dw[0]), 'expected a linear wavelength grid'

        #- top redrock fits (display lists) and models for alternate fits
        fitlists = None
        if zfits is not None:
            #- sky fibers: no alternate fits (redrock fits to sky residuals are noise)
            fitlists = [None if (d is None or str(fibermap_value(sp.fibermap, i, 'OBJTYPE', default='')).strip() == 'SKY')
                        else _fit_table(d)
                        for i, d in enumerate(zfits)]
        alt_model, owner, alt_settings = None, None, None
        if fitlists is not None and alt_models > 0:
            maxfits = max((len(d) for d in zfits if d is not None), default=1)
            nalt = min(alt_models, maxfits - 1)
            alt_model, owner = alt_fit_models(sp, wave, zcat, zfits, fitlists, nalt=nalt, max_dchi2=alt_max_dchi2)
            n = 0 if alt_model is None else len(alt_model)
            log.info(f'models for {n} alternate fits (fits 2-{nalt + 1} with Δχ² < {alt_max_dchi2:g}, not sky/unassigned)')
            if nalt > 0:
                #- the settings, even if no fit qualified, so stack can tell them from "none stored";
                #- JSON has no Infinity: null = no limit
                alt_settings = dict(n=nalt, max_dchi2=float(alt_max_dchi2) if np.isfinite(alt_max_dchi2) else None)

        #- per-pixel arrays, quantized
        arrays = {}
        for name, q in quantize(flux, noise, model, model_precision=model_precision).items():
            arrays[name] = dict(ARRAY_INFO[name], **q)
        if alt_model is not None:
            alt_rows = owner[:, 0]
            arrays['amodel'] = dict(ARRAY_INFO['amodel'], **quantize_model(alt_model, noise[alt_rows], alt_model_precision),
                                    owner=owner)

        #- per-target tables: a subset of the fibermap (DESI names), grouping values, redshifts
        fm = sp.fibermap
        cols = [Column(np.asarray(fm['TARGETID'], dtype=np.int64), name='TARGETID')]
        for name, dtype in GROUP_COLUMNS.items():
            c = _column_from_values(name, groups[name], dtype)
            if c is not None:
                cols.append(c)
        for name in FIBERMAP_COLUMNS + tuple(c for c in fm.colnames if TARGET_BIT_COLUMN.match(c)):
            if name in fm.colnames:
                cols.append(_clean(fm[name]))
        fibermap = Table(cols)

        redshifts = None
        if zcat is not None:
            zcols = [Column(np.asarray(zcat['TARGETID'], dtype=np.int64), name='TARGETID')]
            for name in REDSHIFT_COLUMNS:
                if name in zcat.colnames:
                    c = _clean(zcat[name])
                    if name == 'SUBTYPE':      #- redrock SUBTYPE can be masked ('--')
                        c = Column(np.char.strip(np.asarray(c.filled('') if hasattr(c, 'filled') else c, dtype=str), '-'), name=name)
                    zcols.append(c)
            redshifts = Table(zcols)

        masked_cams = [[[float(sp.wave[b][0]), float(sp.wave[b][-1]), b] for b in bands
                        if not np.any((sp.ivar[b][i] > 0) & ((sp.mask[b][i] == 0) if sp.mask is not None else True))]
                       for i in range(nspec)]
        extra = Table([_objcol('snr', [{b: (round(float(snr[b][i]), 2) if np.isfinite(snr[b][i]) else None)
                                         for b in sp.bands} for i in range(nspec)]),
                       _objcol('masked_cams', masked_cams)])
        if fitlists is not None:
            extra.add_column(_objcol('zfits', fitlists))
        extra['src'] = np.zeros(nspec, dtype=np.int32)

        sources = [dict(source=provenance.get('source'),
                        redrock=provenance.get('redrock') if zcat is not None else None,
                        model_source=model_source,
                        rrdetails=provenance.get('rrdetails'))]
        return cls(fibermap, redshifts, extra, arrays, wave[0], dw[0], len(wave), sources,
                   spgrp=spgrp, overlaps=overlaps, alt_models=alt_settings, title=title)

    @classmethod
    def from_payload(cls, d):
        """
        PageData from a page's payload dict.  Raises ValueError for a format
        version this code can't read, or for keys it doesn't know.
        """
        if check_readable(d.get('format_version')):
            #- a newer minor: keep what this code knows, which is all of format FORMAT_VERSION
            d, dropped = _drop_newer(d)
            msg = f'page data format {d["format_version"]} is newer than this prospect2 reads ({FORMAT_VERSION})'
            if dropped:
                msg += f'; dropping what it doesn\'t know ({", ".join(dropped)}); upgrade prospect2 to keep them'
            get_logger().warning(msg)
        #- Unknown top-level and lowercase per-target keys are ignored: up to FORMAT_VERSION they can
        #- only be derived values (added without a version bump, recomputed when written).  UPPERCASE
        #- keys are stored columns and must be in `columns`.
        n, nwave = d['nspec'], d['nwave']
        columns, targets = d['columns'], d['targets']
        stored = set(columns['fibermap']) | set(columns.get('redshifts') or {})
        unknown = sorted(k for k in targets if k == k.upper() and k not in stored)
        if unknown:
            raise ValueError(f'per-target columns {unknown} are not in the payload\'s columns map (format {d["format_version"]})')
        for k in stored | set(STORED_KEYS):        #- only the keys read here (others are ignored)
            if k in targets and len(targets[k]) != n:
                raise ValueError(f'targets[{k!r}] has {len(targets[k])} entries for {n} spectra')
        fibermap = Table([_column_from_json(c, targets[c], code) for c, code in columns['fibermap'].items()])
        redshifts = None
        if columns.get('redshifts') is not None:
            redshifts = Table([Column(np.array(fibermap['TARGETID']), name='TARGETID')]
                              + [_column_from_json(c, targets[c], code) for c, code in columns['redshifts'].items()])
        extra = Table([_objcol(k, targets[k]) for k in STORED_KEYS if k in targets and k != 'src'])
        extra['src'] = np.array(targets.get('src', [0] * n), dtype=np.int32)
        arrays = {}
        for name, a in d['arrays'].items():
            if a.get('pack') != 'esc8':
                raise ValueError(f'array {name}: unknown packing {a.get("pack")!r}')
            i8, esc = unpack(a['data'], nwave)
            arr = dict(rows=a['rows'], enc=a['enc'], i8=i8, esc=esc)
            for k, v in a.items():
                if k not in ('rows', 'pack', 'enc', 'data'):
                    if len(v) != len(i8):
                        raise ValueError(f'array {name}: {k} has {len(v)} entries for {len(i8)} rows')
                    arr[k] = np.array(v, dtype=np.int64 if k == 'owner' else None)
            arrays[name] = arr
        return cls(fibermap, redshifts, extra, arrays, d['wave0'], d['dwave'], nwave, d['sources'],
                   spgrp=d['spgrp'], overlaps=d['overlaps'], alt_models=d['alt_models'],
                   flux_unit=d['flux_unit'], title=d['title'])

    @classmethod
    def from_html(cls, html):
        """PageData from the HTML of a prospect2 page (format 1.0 or later)"""
        m = re.search(r'<script type="application/json" id="prospect2-data">', html)
        if m is None:
            raise ValueError('not a prospect2 page of format >= 1.0 (no prospect2-data element); regenerate it')
        payload, _ = json.JSONDecoder().raw_decode(html, m.end())
        return cls.from_payload(payload)

    # ------------------------------------------------------------ output

    def _check_rows(self):
        """ValueError if fibermap or redshifts were reordered in place (their rows must match the arrays)"""
        for name, tab in (('fibermap', self.fibermap), ('redshifts', self.redshifts)):
            if tab is not None and (len(tab) != len(self._targetids) or
                                    np.any(np.asarray(tab['TARGETID'], dtype=np.int64) != self._targetids)):
                raise ValueError(f'{name} rows no longer match the spectra (sorted or edited in place?); '
                                 'reorder a page with p[order], e.g. p[np.argsort(p.redshifts["Z"])]')

    def _groups(self):
        """grouping values per row (grouping.grouping_columns format) from the fibermap"""
        n = len(self)
        out = dict(SPGRP=[self.spgrp] * n)
        for name in GROUP_COLUMNS:
            if name in self.fibermap.colnames:
                out[name] = [fibermap_value(self.fibermap, i, name) for i in range(n)]
                out[name] = [None if v is None else (str(v) if GROUP_COLUMNS[name] is str else int(v)) for v in out[name]]
            else:
                out[name] = [None] * n
        return out

    def _derived(self):
        """per-target display values, recomputed from the tables"""
        from redrock.zwarning import ZWarningMask
        try:
            from desispec.maskbits import fibermask
        except ImportError:
            fibermask = None
        fm, n = self.fibermap, len(self)
        col = lambda i, *names, default=None: fibermap_value(fm, i, *names, default=default)
        out = {k: [] for k in DERIVED_KEYS if k != 'zwarn_names'}
        for i in range(n):
            classes, primary = target_classes(fm, i)
            fiberstatus = int(col(i, 'COADD_FIBERSTATUS', 'FIBERSTATUS', default=0))
            out['classes'].append(classes)
            out['primary'].append(primary)
            out['fiberstatus'].append(fiberstatus)
            out['fiberstatus_names'].append(fibermask.names(fiberstatus) if fibermask else [])
            out['exptime'].append(_num(col(i, 'COADD_EXPTIME', 'EXPTIME'), 1))
            out['numexp'].append(int(col(i, 'COADD_NUMEXP', default=1)))
        if self.redshifts is not None and 'ZWARN' in self.redshifts.colnames:
            names = [k for k in dir(ZWarningMask) if k.isupper()]
            out['zwarn_names'] = [[k for k in names if zw & getattr(ZWarningMask, k)]
                                  for zw in (int(fibermap_value(self.redshifts, i, 'ZWARN', default=0)) for i in range(n))]
        return out

    def payload(self, title=None):
        """
        The page's data as a JSON-serializable dict (keys in doc/implementation.md
        "Payload contents").  `title` default: self.title, else the single
        source's file name, else 'prospect2'.
        """
        n = len(self)
        self._check_rows()
        names = [s.get('source') for s in self.sources]
        single = names[0] if len(names) == 1 else None
        title = title or self.title or single or 'prospect2'
        spgrp, keycols = key_columns(self._groups(), self.fibermap['TARGETID'])
        #- VI keys: the input file for a single-source page (so VI survives regenerating,
        #- selecting or restacking it), else the title.  ROW keys only mean something for this
        #- exact list of spectra, so then also a hash of its spectrum identities (TARGETID + grouping
        #- values: TARGETID alone doesn't tell two observations of one target apart).
        vi_key = single or title
        if keycols == ['ROW']:
            ids = ';'.join(','.join('' if v is None else str(v) for v in ident) for ident in _identities(self))
            vi_key += ':rows-' + hashlib.sha1(ids.encode()).hexdigest()[:10]

        targets, columns = {}, dict(fibermap={}, redshifts=None)
        for tname, tab in (('fibermap', self.fibermap), ('redshifts', self.redshifts)):
            if tab is None:
                continue
            columns[tname] = {}
            for c in tab.colnames:
                if tname == 'redshifts' and c == 'TARGETID':
                    continue
                if c != c.upper() or c in targets:
                    raise ValueError(f'{tname} column {c!r}: page columns must be UPPERCASE and unique across fibermap and redshifts')
                targets[c], columns[tname][c] = _json_column(tab[c])
        for k in STORED_KEYS:
            if k in self.extra.colnames:
                targets[k] = [_plain(v) for v in self.extra[k]]
        targets.update(self._derived())

        arrays = {}
        for name, a in self.arrays.items():
            out = dict(rows=a['rows'], pack='esc8', enc=a['enc'], data=pack(a['i8'], a['esc']))
            for k in _side_keys(a):
                out[k] = _plain(np.asarray(a[k]))
            arrays[name] = out

        return dict(
            #- always the current version: the payload is written by this code (derived keys
            #- recomputed; stored keys an older input lacks are tolerated by every reader and
            #- template of this major)
            format_version=FORMAT_VERSION,
            generator=f'prospect2 {__version__}',
            title=title,
            vi_key=vi_key,
            sources=[{k: s.get(k) for k in SOURCE_KEYS} for s in self.sources],
            spgrp=spgrp, keycols=keycols,
            wave0=self.wave0, dwave=self.dwave, nwave=self.nwave, nspec=n,
            flux_unit=self.flux_unit,
            columns=columns,
            targets=targets,
            lines=read_lines(),
            overlaps=self.overlaps,
            alt_models=self.alt_models,
            arrays=arrays,
        )

    def to_html(self, title=None):
        """HTML of the self-contained viewer page, with the current template"""
        from .page import render_html
        if len(self) == 0:
            raise ValueError('no spectra: a page needs at least one')
        return render_html(self.payload(title))


def _owner(v):
    """
    An array's `owner` side list as int64 [nrow, ncol].  Column 0 is always the
    spectrum the row belongs to (slicing and stacking use only that column); the
    other columns' meaning depends on the row space (alt: index into zfits).
    """
    v = np.asarray(v, dtype=np.int64)
    if v.size == 0:
        return v.reshape(0, 2)
    if v.ndim != 2:
        raise ValueError(f'owner must be [nrow, ncol], got shape {v.shape}')
    return v


def _plain(v):
    """numpy scalars/arrays inside a value -> plain Python, for JSON"""
    if isinstance(v, np.ndarray):
        return [_plain(x) for x in v]
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, dict):
        return {k: _plain(x) for k, x in v.items()}
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return float(v) if np.isfinite(v) else None       #- JSON has no NaN/Infinity
    return v


def stack(pages, title=None):
    """
    Stack PageData objects into one, like desispec.spectra.stack: tables are
    vstacked (a column missing from some pages is masked there), per-pixel
    arrays concatenated losslessly, sources combined.  `title` is the page
    title (the label at the top left), stored as the result's .title; by
    default the pages' titles are not kept, and a written page is titled
    after its single source file, else write_page's file name.

    Raises ValueError if the pages differ in wavelength grid, model
    availability, alternate-fit settings or spectral grouping, or if the stack
    would hold the same spectrum twice: TARGETID and every known grouping
    value equal (a value unknown on one page, e.g. NIGHT, counts as different;
    if that makes the spectra not unique by their key, VI falls back to rows,
    with a warning when written).  Empty pages (e.g. a selection that
    matched nothing) are fine; writing an empty page is not.
    """
    pages = list(pages)
    if not pages:
        raise ValueError('no pages to stack')
    for p in pages:
        p._check_rows()
    p0 = pages[0]
    spgrps = sorted({p.spgrp for p in pages if p.spgrp is not None})
    if len(spgrps) > 1:
        raise ValueError(f'pages mix spectral groupings (SPGRP): {", ".join(spgrps)}')
    #- optional parts (top redrock fits, alternate-fit models): a page without them adds no rows
    optional = {'amodel'}
    #- alternate-fit model settings matter only for pages that show alternate fits (zfits);
    #- among those, None (no models stored) is a setting of its own
    has_fits = lambda p: 'zfits' in p.extra.colnames and any(z is not None for z in p.extra['zfits'])
    alt_settings = [p.alt_models for p in pages if has_fits(p)]
    for j, p in enumerate(pages[1:], 1):
        for what, a, b in (('wavelength grid', (p0.wave0, p0.dwave, p0.nwave), (p.wave0, p.dwave, p.nwave)),
                           ('flux unit', p0.flux_unit, p.flux_unit),
                           ('camera overlaps', p0.overlaps, p.overlaps),
                           ('arrays (e.g. some pages have no model)', sorted(set(p0.arrays) - optional),
                            sorted(set(p.arrays) - optional)),
                           ('redshifts (some pages have none)', p0.redshifts is None, p.redshifts is None),
                           ('alternate-fit model settings', alt_settings[0] if alt_settings else None,
                            p.alt_models if has_fits(p) else (alt_settings[0] if alt_settings else None))):
            if a != b:
                raise ValueError(f'can not stack: page {j} differs from page 0 in {what}: {b} vs {a}')
    names = [n for n in dict.fromkeys(n for p in pages for n in p.arrays)]
    extras = [c for c in dict.fromkeys(c for p in pages for c in p.extra.colnames)]
    template = {n: next(p.arrays[n] for p in pages if n in p.arrays) for n in names}
    pages = [_with_optional(p, names, extras, template) for p in pages]
    p0 = pages[0]

    #- sources: combined, identical ones merged
    sources, src = [], []
    for p in pages:
        remap = []
        for s in p.sources:
            if s not in sources:
                sources.append(s)
            remap.append(sources.index(s))
        src.append(np.array([remap[int(i)] for i in p.extra['src']], dtype=np.int32))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')        #- metadata merge warnings
        fibermap = vstack([p.fibermap for p in pages], join_type='outer', metadata_conflicts='silent')
        redshifts = None if p0.redshifts is None else vstack([p.redshifts for p in pages], join_type='outer',
                                                              metadata_conflicts='silent')
        extra = vstack([p.extra for p in pages], metadata_conflicts='silent')
    extra['src'] = np.concatenate(src)

    arrays, offset = {}, np.cumsum([0] + [len(p) for p in pages])
    for name, a0 in p0.arrays.items():
        parts = [p.arrays[name] for p in pages]
        a = dict(rows=a0['rows'], enc=a0['enc'],
                 i8=np.concatenate([x['i8'] for x in parts]), esc=np.concatenate([x['esc'] for x in parts]).astype('<i4'))
        for k in _side_keys(a0):
            vals = [np.asarray(x[k]) for x in parts]
            if k == 'owner':
                vals = [_owner(v).copy() for v in vals]      #- copy: don't modify the inputs
                for v, off in zip(vals, offset):
                    v[:, 0] += off
            a[k] = np.concatenate(vals)
        arrays[name] = a

    out = PageData(fibermap, redshifts, extra, arrays, p0.wave0, p0.dwave, p0.nwave, sources,
                   spgrp=spgrps[0] if spgrps else None, overlaps=p0.overlaps,
                   alt_models=alt_settings[0] if alt_settings else None,
                   flux_unit=p0.flux_unit, title=title)

    #- the same spectrum twice (e.g. one page stacked with itself) is an error, unless the
    #- inputs already had those duplicates
    ndup = lambda ids: len(ids) - len(set(ids))
    if ndup(_identities(out)) > sum(ndup(_identities(p)) for p in pages):
        raise ValueError('stacked pages contain the same spectrum more than once')
    return out


def _identities(p):
    """per spectrum: (TARGETID + all grouping values), which identifies a spectrum"""
    g = p._groups()
    return list(zip([int(t) for t in p.fibermap['TARGETID']], *(g[c] for c in GROUP_COLUMNS)))


def _with_optional(p, names, extras, template):
    """
    Shallow copy of page `p` with zero-row arrays for optional arrays it lacks
    (shaped like `template[name]`) and None-filled extra columns, so it stacks
    with pages that have them.
    """
    if all(n in p.arrays for n in names) and all(c in p.extra.colnames for c in extras):
        return p
    arrays = dict(p.arrays)
    for n in names:
        if n not in arrays:
            t = template[n]
            arrays[n] = dict(rows=t['rows'], enc=t['enc'], i8=np.zeros((0, p.nwave), dtype='i1'),
                             esc=np.zeros(0, dtype='<i4'))
            for k in _side_keys(t):
                arrays[n][k] = np.zeros((0,) + np.asarray(t[k]).shape[1:], dtype=np.asarray(t[k]).dtype)
    extra = p.extra.copy(copy_data=False)
    for c in extras:
        if c not in extra.colnames:
            extra.add_column(_objcol(c, [None] * len(p)), index=len(extra.colnames) - 1)
    return PageData(p.fibermap, p.redshifts, extra, arrays, p.wave0, p.dwave, p.nwave, p.sources,
                    spgrp=p.spgrp, overlaps=p.overlaps, alt_models=p.alt_models, flux_unit=p.flux_unit,
                    title=p.title)

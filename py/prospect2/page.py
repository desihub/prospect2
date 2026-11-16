#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.page
==============

Build a prospect2 viewer page from in-memory objects.  The main entry point is
make_html(), which takes desispec Spectra (and optionally redrock outputs) and
returns the HTML of one self-contained page.
"""

import os, re, json, html, warnings
import numpy as np
from astropy.table import vstack

from desiutil.log import get_logger
from desispec.coaddition import coadd_cameras

from .encode import encode, encode_model
from .models import camera_models, merge_camera_models, alt_fit_models
from .io import DATADIR, read_lines
from .util import fibermap_value
from .grouping import grouping_columns, concat_columns, key_columns

#- redrock's "same redshift" criterion for marking near-duplicate alternate fits
_DUP_KMS = 1000


def mag(flux):
    """AB magnitude of a nanomaggie flux, or None if flux is None or not > 0"""
    if flux is None:
        return None
    flux = float(flux)
    return round(22.5 - 2.5 * np.log10(flux), 2) if flux > 0 else None


def _num(x, ndigits):
    """round(float(x), ndigits), or None if x is None or not finite"""
    if x is None or not np.isfinite(float(x)):
        return None
    return round(float(x), ndigits)


def target_classes(row):
    try:
        from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
    except ImportError:
        return [], []
    names = []
    for col, mask in [('DESI_TARGET', desi_mask), ('BGS_TARGET', bgs_mask),
                      ('MWS_TARGET', mws_mask), ('SCND_TARGET', scnd_mask)]:
        if col in row.colnames and row[col] != 0:
            names += [n for n in mask.names(int(row[col]))
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


def _combine(spectra, zcat, models, fits):
    """
    Stack lists of inputs (one entry per input file) into single objects,
    raising ValueError if the redrock template versions differ (models computed
    from templates must use the templates of the fit).  Also returns the
    per-row grouping values (grouping.grouping_columns), read before stacking
    because stack keeps only the first input's headers.
    """
    from desispec.spectra import stack
    n = len(spectra)
    for name, x in (('zcat', zcat), ('models', models), ('fits', fits)):
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
    sp = stack(list(spectra))
    if zcat is not None:
        zcat = vstack(list(zcat), metadata_conflicts='silent')
    if models is not None:
        models = {b: np.concatenate([m[b] for m in models]) for b in models[0]}
    if fits is not None:
        fits = [f for flist in fits for f in flist]
    return sp, zcat, models, fits, groups


def make_payload(spectra, zcat=None, models=None, fits=None, title=None, provenance=None,
                 alt_models=2, alt_max_dchi2=100, model_precision=100, alt_model_precision=30,
                 groups=None):
    """
    The page's data payload (a JSON-serializable dict); see make_html for the
    arguments.  Inputs must already be single objects (not lists).  `groups`
    (internal) are the per-row grouping values of stacked inputs; default
    grouping.grouping_columns(spectra).  The keys,
    and how to add one, are listed in doc/implementation.md "Payload contents";
    keep that in sync.
    """
    log = get_logger()
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
    spgrp, keycols = key_columns(groups, tids)
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
    if fits is not None:
        if zcat is None:
            raise ValueError('fits given without zcat')
        if len(fits) != nspec:
            raise ValueError(f'fits has {len(fits)} entries for {nspec} spectra')

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

    #- Per-target metadata
    from redrock.zwarning import ZWarningMask
    zw_names = [n for n in dir(ZWarningMask) if n.isupper()]
    try:
        from desispec.maskbits import fibermask
    except ImportError:
        fibermask = None

    fitlists = None
    if fits is not None:
        #- sky fibers: no alternate fits (redrock fits to sky residuals are noise)
        fitlists = [None if (d is None or str(fibermap_value(sp.fibermap, i, 'OBJTYPE', default='')).strip() == 'SKY')
                    else _fit_table(d)
                    for i, d in enumerate(fits)]

    alt_model, alt_rows = None, None
    if fitlists is not None and alt_models > 0:
        maxfits = max((len(d) for d in fits if d is not None), default=1)
        nalt = min(alt_models, maxfits - 1)
        alt_model, alt_rows = alt_fit_models(sp, wave, zcat, fits, fitlists, nalt=nalt, max_dchi2=alt_max_dchi2)
        n = 0 if alt_model is None else len(alt_model)
        log.info(f'models for {n} alternate fits (fits 2-{nalt + 1} with Δχ² < {alt_max_dchi2:g}, not sky/unassigned)')

    #- only TARGETID is required: other columns fall back (to None, '', or the
    #- per-exposure EXPTIME/FIBERSTATUS and NUMEXP=1 for uncoadded spectra)
    fm = sp.fibermap
    targets = []
    for i in range(nspec):
        r = fm[i]
        col = lambda *names, default=None: fibermap_value(fm, i, *names, default=default)
        allclass, primary = target_classes(r)
        fiberstatus = int(col('COADD_FIBERSTATUS', 'FIBERSTATUS', default=0))
        t = dict(
            i=i, targetid=str(int(r['TARGETID'])),
            fiber=groups['FIBER'][i], tileid=groups['TILEID'][i], night=groups['NIGHT'][i],
            expid=groups['EXPID'][i], survey=groups['SURVEY'][i] or '', program=groups['PROGRAM'][i] or '',
            ra=_num(col('TARGET_RA'), 6), dec=_num(col('TARGET_DEC'), 6),
            objtype=str(col('OBJTYPE', default='')), morph=str(col('MORPHTYPE', default='')).strip(),
            classes=allclass, primary=primary,
            mag={k: mag(col(f'FLUX_{k.upper()}')) for k in ('g', 'r', 'z', 'w1', 'w2')},
            fibermag_r=mag(col('FIBERFLUX_R')),
            ebv=_num(col('EBV'), 3),
            exptime=_num(col('COADD_EXPTIME', 'EXPTIME'), 1), numexp=int(col('COADD_NUMEXP', default=1)),
            fiberstatus=fiberstatus,
            fiberstatus_names=(fibermask.names(fiberstatus) if fibermask else []),
            snr={b: (round(float(snr[b][i]), 2) if np.isfinite(snr[b][i]) else None) for b in sp.bands},
            #- cameras with no unmasked pixels for this spectrum: [lo, hi, band]
            masked_cams=[[float(sp.wave[b][0]), float(sp.wave[b][-1]), b] for b in bands
                         if not np.any((sp.ivar[b][i] > 0) & ((sp.mask[b][i] == 0) if sp.mask is not None else True))],
        )
        if zcat is not None:
            zb = zcat[i]
            zwarn = int(zb['ZWARN'])
            t.update(z=float(zb['Z']), zerr=float(zb['ZERR']), zwarn=zwarn,
                     zwarn_names=[n for n in zw_names if zwarn & getattr(ZWarningMask, n)],
                     spectype=str(zb['SPECTYPE']).strip(), subtype=str(zb['SUBTYPE']).strip().strip('-'),
                     deltachi2=round(float(zb['DELTACHI2']), 2), chi2=round(float(zb['CHI2']), 1),
                     npixels=int(zb['NPIXELS']))
        if fitlists is not None and fitlists[i] is not None:
            t['fits'] = fitlists[i]
        targets.append(t)

    #- the JS keys stored VI and names the VI export by `source`, so it must not be null
    title = title or provenance.get('source') or 'prospect2'
    source = provenance.get('source') or title
    payload = dict(
        title=title,
        source=source,
        redrock=provenance.get('redrock') if zcat is not None else None,
        model_source=model_source,
        rrdetails=provenance.get('rrdetails'),
        spgrp=spgrp, keycols=keycols,
        wave0=float(wave[0]), dwave=float(dw[0]), nwave=len(wave), nspec=nspec,
        flux_unit='1e-17 erg/s/cm²/Å',
        targets=targets,
        lines=read_lines(),
        overlaps=overlaps,
        **encode(flux, noise, model, model_precision=model_precision),
    )
    if alt_model is not None:
        payload['amodel'], payload['amscale'] = encode_model(alt_model, noise[alt_rows], alt_model_precision)
        payload['alt_models'] = dict(n=nalt, max_dchi2=alt_max_dchi2)
    return payload


def render_html(payload):
    """
    Fill the viewer template with uPlot and `payload`; returns the page HTML.
    The title is HTML-escaped, and every "<" in the payload JSON is written as
    \\u003c (same string in JS), so no string can end the <script> element
    ("</script>") or open a comment ("<!--").  The data goes in last, so text in
    it is never mistaken for a placeholder.
    """
    template = (DATADIR / 'prospect2_template.html').read_text(encoding='utf-8')
    uplot_js = (DATADIR / 'vendor' / 'uPlot.iife.min.js').read_text(encoding='utf-8')
    uplot_css = (DATADIR / 'vendor' / 'uPlot.min.css').read_text(encoding='utf-8')
    data = json.dumps(payload, separators=(',', ':')).replace('<', '\\u003c')
    return (template.replace('__TITLE__', html.escape(payload['title']))
                    .replace('/*__UPLOT_CSS__*/', uplot_css)
                    .replace('/*__UPLOT_JS__*/', uplot_js)
                    .replace('"__DATA__"', data))


def make_html(spectra, zcat=None, models=None, fits=None, title=None, provenance=None,
              alt_models=2, alt_max_dchi2=100, model_precision=100, alt_model_precision=30):
    """
    HTML of a self-contained prospect2 viewer page.

    Args:
        spectra: desispec Spectra with per-camera data (not merged across
            cameras), e.g. from read_spectra(); or a list of them, which are
            stacked after checking that their redrock template versions agree.
            When passing one Spectra built by stacking several inputs yourself,
            you are responsible for that consistency, and only the first
            input's headers (SPGRP, NIGHT, SURVEY, PROGRAM) survive the
            stack, so pass the list instead.  Spectra are keyed for VI by
            TILEID/FIBER/NIGHT(/EXPID) for tile-based groupings (SPGRP
            cumulative, pernight, perexp) and by SURVEY/PROGRAM/TARGETID for
            healpix/uniqpix; mixing groupings raises ValueError.

    Options (each row-matched to spectra; for a list of spectra, a list of the same length):
        zcat: redrock REDSHIFTS table; default spectra.redshifts
            (read_spectra(..., return_redshifts=True)).  Its meta should hold the
            TEMNAMnn/TEMVERnn keywords so computed models use the fit's templates.
        models: per-camera best-fit models {band: [nspec, nwave_band]}, e.g. from
            prospect2.io.read_rrmodel; default spectra.model
            (read_spectra(..., return_models=True)); if neither, computed from
            the redrock templates + COEFF.
        fits: alternate redrock fits from prospect2.io.read_rrdetails.
        title: page title; default provenance['source'].
        provenance: file names shown on the page, keys source, redrock,
            model_source, rrdetails.
        alt_models: store models for alternate fits 2..alt_models+1 (0 = none).
        alt_max_dchi2: only for alternate fits with chi2 - best chi2 below this.
        model_precision, alt_model_precision: store models in steps of
            (median data noise)/N.

    Explicit zcat/models override spectra.redshifts/spectra.model.  Raises
    ValueError if the inputs are not row-matched or consistent.
    """
    if isinstance(spectra, (list, tuple)):
        spectra, zcat, models, fits, groups = _combine(spectra, zcat, models, fits)
    else:
        groups = None
    payload = make_payload(spectra, zcat=zcat, models=models, fits=fits, title=title, groups=groups,
                           provenance=provenance, alt_models=alt_models, alt_max_dchi2=alt_max_dchi2,
                           model_precision=model_precision, alt_model_precision=alt_model_precision)
    return render_html(payload)


def write_html(html, outfile):
    """Write page `html` to `outfile`, logging whether it replaced an existing file"""
    log = get_logger()
    verb = 'Overwrote existing' if os.path.exists(outfile) else 'Wrote'
    with open(outfile, 'w', encoding='utf-8') as fx:
        fx.write(html)
    log.info(f'{verb} {outfile}: {os.path.getsize(outfile)/1e6:.2f} MB')

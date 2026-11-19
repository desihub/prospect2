#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.io
============

Finding and reading the redrock files that go with a coadd, and the package data.
"""

import os
from importlib.resources import files

import numpy as np
from astropy.table import Table
from astropy.io import fits

from desiutil.log import get_logger
from desispec.io.util import checkgzip

#- package data: viewer template, vendored uPlot, line lists
DATADIR = files('prospect2') / 'data'


def find_file(filename):
    """filename or filename(.gz) if it exists, else None"""
    try:
        return checkgzip(filename)
    except (FileNotFoundError, IOError):
        return None


def find_redrock_files(infile):
    """
    Redrock files matching a coadd or spectra file, as a dict with keys
    redrock, rrmodel, rrdetails: coadd-X.fits -> redrock-X.fits,
    rrmodel-X.fits, rrdetails-X.h5 (.gz optional).  Missing files are None,
    and rrmodel/rrdetails are only looked for next to a redrock file.
    """
    basedir, basefile = os.path.split(infile)
    rrfile = basefile.replace('coadd', 'redrock').replace('spectra', 'redrock')
    out = dict(redrock=find_file(os.path.join(basedir, rrfile)), rrmodel=None, rrdetails=None)
    if out['redrock']:
        rrdir, rrbase = os.path.split(out['redrock'])
        if rrbase.startswith('redrock'):
            suffix = rrbase[len('redrock'):]
            out['rrmodel'] = find_file(os.path.join(rrdir, 'rrmodel' + suffix))
            out['rrdetails'] = find_file(os.path.join(rrdir, 'rrdetails' + suffix.replace('.fits.gz', '.fits').replace('.fits', '.h5')))
    return out


def read_redshifts(rrfile):
    """
    Redrock REDSHIFTS table, with the primary header (incl. the TEMNAMnn/TEMVERnn
    template versions) merged into its meta.
    """
    zcat = Table.read(rrfile, 'REDSHIFTS')
    zcat.meta.update(fits.getheader(rrfile, 0))
    return zcat


def read_rrmodel(rrmodel_file, targetids, wave):
    """
    Read per-camera redrock models {band: [nspec, nwave_band]} from an rrmodel
    file (REDSHIFTS + {B,R,Z}_WAVELENGTH + {B,R,Z}_MODEL HDUs), for `targetids`
    in that order.  `wave` is the spectra's {band: wavelength}.  Raises ValueError
    if the file doesn't cover every TARGETID and camera on the same wavelength grids.
    """
    with fits.open(rrmodel_file, memmap=True) as hx:
        names = {h.name for h in hx}
        tids = np.asarray(hx['REDSHIFTS'].data['TARGETID'])
        row = {t: i for i, t in enumerate(tids)}
        missing = [t for t in targetids if t not in row]
        if missing:
            raise ValueError(f'{len(missing)} TARGETIDs (e.g. {missing[0]}) not in {rrmodel_file}')
        idx = np.array([row[t] for t in targetids])
        out = {}
        for b in wave:
            B = b.upper()
            if f'{B}_MODEL' not in names or f'{B}_WAVELENGTH' not in names:
                raise ValueError(f'no {B}_MODEL/{B}_WAVELENGTH in {rrmodel_file}')
            w = hx[f'{B}_WAVELENGTH'].data
            if len(w) != len(wave[b]) or not np.allclose(w, wave[b]):
                raise ValueError(f'{B}_WAVELENGTH in {rrmodel_file} differs from the spectra')
            out[b] = hx[f'{B}_MODEL'].data[idx].astype(np.float64)
    return out


def read_rrdetails(rrdetails_file, zcat, nfits=5):
    """
    Top-`nfits` redrock fits per target from an rrdetails*.h5 file: a list with
    one entry per zcat row, each the raw zfit records (numpy structured array
    with z, spectype, subtype, chi2, zwarn, coeff, ...) ranked by chi2, or None
    for a target that is missing, has no fit, or whose best fit doesn't match
    zcat (z and SPECTYPE).  Lists from several files can be concatenated, like
    their zcats.
    """
    import h5py
    log = get_logger()
    s = lambda x: (x.decode() if isinstance(x, bytes) else str(x)).strip()
    out = []
    nbad = 0
    with h5py.File(rrdetails_file, 'r') as f:
        zf = f['zfit']
        for row in zcat:
            key = str(int(row['TARGETID']))
            if key not in zf or row['NPIXELS'] == 0:
                out.append(None)
                continue
            d = zf[key]['zfit'][:nfits]
            if abs(d['z'][0] - row['Z']) > 1e-6 or s(d['spectype'][0]) != s(row['SPECTYPE']):
                nbad += 1
                out.append(None)
                continue
            out.append(d)
    if nbad:
        log.warning(f'{nbad} targets in {rrdetails_file} have a best fit that does not match the redrock file; alternate fits not shown for them')
    return out


def air_to_vac(air):
    s2 = (1e4 / air) ** 2
    n = 1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - s2) \
        + 1.599740894897e-4 / (3.892568793293e1 - s2)
    return air * n


def read_lines():
    lines = []
    for kind in ('emission', 'absorption'):
        fn = DATADIR / f'{kind}_lines.txt'
        for row in fn.read_text(encoding='utf-8').splitlines():
            if row.startswith('#') or not row.strip():
                continue
            name, longname, lam, vac, major = [x.strip() for x in row.split(',')]
            lam = float(lam)
            if vac != 'True':
                lam = air_to_vac(lam)
            lines.append(dict(name=name, longname=longname, wave=round(lam, 2),
                              kind=kind[:3], major=(major == 'True')))
    return lines


def write_html(html, outfile):
    """Write page `html` to `outfile`, logging whether it replaced an existing file"""
    log = get_logger()
    verb = 'Overwrote existing' if os.path.exists(outfile) else 'Wrote'
    with open(outfile, 'w', encoding='utf-8') as fx:
        fx.write(html)
    log.info(f'{verb} {outfile}: {os.path.getsize(outfile)/1e6:.2f} MB')


def read_page(filename):
    """
    PageData from a prospect2 page (cf. desispec.io.read_spectra).  Raises
    ValueError if the page's data format can't be read by this code.
    """
    from .pagedata import PageData
    with open(filename, encoding='utf-8') as fx:
        return PageData.from_html(fx.read())


def write_page(filename, page, title=None):
    """
    Write PageData `page` as a viewer page (cf. desispec.io.write_spectra),
    creating the directory if needed.  `title` default: page.title, else its
    single source file, else (stacked pages) the output file's base name.
    """
    log = get_logger()
    if title is None and page.title is None and len(page.sources) > 1:
        title = os.path.splitext(os.path.basename(filename))[0]
    html = page.to_html(title=title)          #- before creating directories: it may raise
    outdir = os.path.dirname(filename)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)
        log.info(f'Created directory {outdir}')
    write_html(html, filename)

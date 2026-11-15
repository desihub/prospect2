#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.test.util
===================

Helpers shared by the tests: synthetic inputs, a Python decoder for the "q2"
payload encoding (mirrors decodeQ2() in the viewer template), extracting the
payload from a page, and locating the optional real test data.
"""

import os, gzip, json, base64
import numpy as np
from astropy.table import Table

#- real test data: petal 4 of tile 1496 (coadd, redrock, and if present rrmodel, rrdetails)
TEST_TILE_FILE = 'coadd-4-1496-thru20240214.fits'
REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

#- TARGETIDs above 2**53 check that they reach the page as exact strings
TARGETID0 = 39628266889283065

#- redrock template versions written into synthetic zcat headers
TEMPLATE_META = {'TEMNAM00': 'GALAXY', 'TEMVER00': '2.6', 'TEMNAM01': 'QSO:::LOZ', 'TEMVER01': '1.0'}


def find_test_data():
    """
    Directory holding the real test data (TEST_TILE_FILE and its redrock files):
    $PROSPECT2_TEST_DATA if set, else the repo's scratch/data/, else the tile
    directory in the loa or nevis production under $DESI_ROOT.  Returns None if
    none of them has the data.  Raises FileNotFoundError if $PROSPECT2_TEST_DATA
    is set but doesn't have it, rather than silently testing something else.
    """
    envdir = os.environ.get('PROSPECT2_TEST_DATA')
    if envdir:
        if not os.path.exists(os.path.join(envdir, TEST_TILE_FILE)):
            raise FileNotFoundError(f'$PROSPECT2_TEST_DATA={envdir} has no {TEST_TILE_FILE}')
        return envdir
    candidates = [os.path.join(REPO_DIR, 'scratch', 'data')]
    if 'DESI_ROOT' in os.environ:
        for prod in ('loa', 'nevis'):
            candidates.append(os.path.join(os.environ['DESI_ROOT'], 'spectro', 'redux', prod,
                                           'tiles', 'cumulative', '1496', '20240214'))
    for d in candidates:
        if os.path.exists(os.path.join(d, TEST_TILE_FILE)):
            return d
    return None


def synthetic_spectra(nspec=4, seed=1):
    """
    desispec Spectra with realistic b/r/z grids, unit-ish flux with noise
    sigma=0.5, no masked pixels, and a fibermap with the columns the page reads.
    TARGETIDs start at TARGETID0, fibers at 0.
    """
    from desispec.test.util import get_blank_spectra
    rng = np.random.default_rng(seed)
    sp = get_blank_spectra(nspec)
    for b in sp.bands:
        shape = sp.flux[b].shape
        sp.ivar[b][:] = 4.0
        sp.flux[b][:] = 1.0 + rng.normal(0, 0.5, size=shape)
    fm = sp.fibermap
    fm['TARGETID'] = TARGETID0 + np.arange(nspec)
    fm['OBJTYPE'] = 'TGT'
    fm['MORPHTYPE'] = 'REX'
    for col in ('FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_R'):
        fm[col] = 10.0
    fm['TILEID'] = 1496
    fm['COADD_EXPTIME'] = 900.0
    fm['COADD_NUMEXP'] = 1
    fm['COADD_FIBERSTATUS'] = 0
    return sp


def synthetic_zcat(sp, z=0.5, meta=None):
    """Redrock-like REDSHIFTS table row-matched to `sp`, with template-version meta"""
    n = sp.num_spectra()
    zcat = Table()
    zcat['TARGETID'] = np.asarray(sp.fibermap['TARGETID'])
    zcat['Z'] = np.full(n, z)
    zcat['ZERR'] = np.full(n, 1e-4)
    zcat['ZWARN'] = np.zeros(n, dtype=np.int64)
    zcat['SPECTYPE'] = np.array(['GALAXY'] * n)
    zcat['SUBTYPE'] = np.array([''] * n)
    zcat['DELTACHI2'] = np.full(n, 100.0)
    zcat['CHI2'] = np.full(n, 8000.0)
    zcat['NPIXELS'] = np.full(n, 7000)
    zcat.meta.update(TEMPLATE_META if meta is None else meta)
    return zcat


def synthetic_models(sp, value=1.0):
    """Per-camera models {band: [nspec, nwave_band]}, constant `value`"""
    return {b: np.full(sp.flux[b].shape, value) for b in sp.bands}


def synthetic_fits(zlist, chi2list, spectype='GALAXY'):
    """One target's raw rrdetails records (structured array), as from read_rrdetails"""
    dt = [('z', 'f8'), ('spectype', 'S6'), ('subtype', 'S20'), ('chi2', 'f8'),
          ('zwarn', 'i8'), ('coeff', 'f8', (10,))]
    d = np.zeros(len(zlist), dtype=dt)
    d['z'] = zlist
    d['chi2'] = chi2list
    d['spectype'] = spectype
    return d


def _unpack(b64, dtype):
    return np.frombuffer(gzip.decompress(base64.b64decode(b64)), dtype=dtype)


def _unesc8(packed, masked_flag=False):
    """Undo _esc8: int8 stream with -128 escapes into an int32 stream (row-major order)"""
    a = _unpack(packed[0], 'i1').astype(np.int64)
    x = _unpack(packed[1], '<i4')
    masked = (a == -127) if masked_flag else None
    out = a.copy()
    esc = (a == -128)
    out[esc] = x
    return out, masked


def decode_q2(payload):
    """
    Python mirror of decodeQ2() in the viewer: returns (flux, noise, model,
    amodel) as float arrays [nspec or nrow, nwave], NaN where masked, None for
    absent models.
    """
    nspec, nw = payload['nspec'], payload['nwave']
    dln, _ = _unesc8(payload['noise'])
    q, masked = _unesc8(payload['flux'], masked_flag=True)
    ln = np.cumsum(dln.reshape(nspec, nw), axis=1)
    noise = np.exp(ln / 100)
    flux = q.reshape(nspec, nw) * noise / 10
    masked = masked.reshape(nspec, nw)
    flux[masked] = np.nan
    noise[masked] = np.nan

    def model(packed, scale):
        if packed is None:
            return None
        dm, _ = _unesc8(packed)
        return np.cumsum(dm.reshape(len(scale), nw), axis=1) * np.asarray(scale)[:, None]

    return flux, noise, model(payload['model'], payload['mscale']), model(payload.get('amodel'), payload.get('amscale'))


def payload_from_html(html):
    """The JSON payload embedded in a generated page"""
    start = html.index('const D = ') + len('const D = ')
    payload, _ = json.JSONDecoder().raw_decode(html, start)
    return payload

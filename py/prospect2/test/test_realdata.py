#
# See top-level LICENSE.rst file for Copyright information

"""
Tier 3: real data (petal 4 of tile 1496), found by prospect2.test.util.find_test_data
($PROSPECT2_TEST_DATA, scratch/data/, or the loa/nevis tile directory under $DESI_ROOT).
Skipped if none exists.  Checks behaviors on specific fibers, not bytes.
"""

import os
import numpy as np
import pytest

from prospect2.test.util import find_test_data, decode_q2, payload_from_html, TEST_TILE_FILE

DATADIR = find_test_data()
pytestmark = pytest.mark.skipif(DATADIR is None, reason='no real test data (set $PROSPECT2_TEST_DATA)')

#- 2000-2019: ordinary targets and skies; 2026-2028: z camera masked (BADAMPZ); 2036: no data;
#- 2183: rrmodel sag next to a masked gap; 2208: camera masked in an overlap; 2277: b+z masked
FIBERS = '2000:2020,2026-2028,2036,2183,2208,2277'
NSPEC = 20 + 3 + 4
#- alternate-fit models need redrock templates
HAVE_TEMPLATES = 'RR_TEMPLATE_DIR' in os.environ


def _file(prefix, ext='.fits'):
    path = os.path.join(DATADIR, TEST_TILE_FILE.replace('coadd', prefix).replace('.fits', ext))
    return path if os.path.exists(path) else None


@pytest.fixture(scope='module')
def page(tmp_path_factory):
    """(payload, html) of a page for FIBERS made by the command-line script"""
    from prospect2.scripts.plot_prospect2 import parse, main
    outfile = str(tmp_path_factory.mktemp('realdata') / 'page.html')
    argv = ['-i', os.path.join(DATADIR, TEST_TILE_FILE), '-f', FIBERS, '-o', outfile]
    if not HAVE_TEMPLATES:
        argv += ['--alt-models', '0']
    assert main(parse(argv)) == 0
    html = open(outfile, encoding='utf-8').read()
    return payload_from_html(html), html


def _target(payload, fiber):
    return next(t for t in payload['targets'] if t['fiber'] == fiber)


def test_masked_fibers(page):
    payload, _ = page
    assert payload['nspec'] == NSPEC
    flux = decode_q2(payload)[0]
    fibers = [t['fiber'] for t in payload['targets']]
    t = _target(payload, 2036)
    assert np.all(np.isnan(flux[fibers.index(2036)]))
    assert t['numexp'] == 0 and [m[2] for m in t['masked_cams']] == ['b', 'r', 'z']
    for fiber in (2026, 2027, 2028):
        assert [m[2] for m in _target(payload, fiber)['masked_cams']] == ['z']
    assert [m[2] for m in _target(payload, 2277)['masked_cams']] == ['b', 'z']
    assert _target(payload, 2208)['masked_cams'] == []


def test_model_source(page):
    payload, _ = page
    assert payload['redrock'] == TEST_TILE_FILE.replace('coadd', 'redrock')
    rrmodel = _file('rrmodel')
    if rrmodel:
        assert payload['model_source'] == os.path.basename(rrmodel)
    else:
        assert payload['model_source'] == 'redrock templates + COEFF'
    #- no model drawn for spectra without a redrock fit (fiber 2036): the JS hides it,
    #- but the stored model must not be NaN anywhere
    assert np.all(np.isfinite(decode_q2(payload)[2]))


@pytest.mark.skipif(_file('rrdetails', '.h5') is None if DATADIR else True, reason='no rrdetails file')
def test_alternate_fits(page):
    payload, _ = page
    with_fits = [t for t in payload['targets'] if 'fits' in t]
    assert with_fits, 'no alternate fits on the page'
    for t in payload['targets']:
        if t['objtype'].strip() == 'SKY':
            assert 'fits' not in t
    if not HAVE_TEMPLATES:
        assert 'amodel' not in payload
        return
    nalt = payload['alt_models']['n']
    nrows = 0
    for t in with_fits:
        for k, f in enumerate(t['fits']):
            if len(f) > 6:                           #- this fit has a stored model
                nrows += 1
                assert 1 <= k <= nalt and f[3] < payload['alt_models']['max_dchi2']
                assert t['objtype'].strip() == 'TGT'
    assert nrows == len(payload['amscale']) > 0


@pytest.mark.skipif(not HAVE_TEMPLATES, reason='needs $RR_TEMPLATE_DIR')
def test_api_matches_script(page):
    """make_html on read_spectra(return_models, return_redshifts) gives the script's page"""
    import prospect2
    from desispec.io import read_spectra
    from prospect2.io import read_rrdetails
    from prospect2.util import select_spectra
    _, html = page
    if _file('rrmodel') is None:
        pytest.skip('needs an rrmodel file')
    sp = read_spectra(os.path.join(DATADIR, TEST_TILE_FILE), return_models=True, return_redshifts=True)
    sp = sp[select_spectra(sp.fibermap, fibers=FIBERS)]
    rrdetails = _file('rrdetails', '.h5')
    fits = read_rrdetails(rrdetails, sp.redshifts) if rrdetails else None
    provenance = dict(source=TEST_TILE_FILE, redrock=TEST_TILE_FILE.replace('coadd', 'redrock'),
                      model_source=os.path.basename(_file('rrmodel')),
                      rrdetails=os.path.basename(rrdetails) if rrdetails else None)
    assert prospect2.make_html(sp, fits=fits, provenance=provenance) == html

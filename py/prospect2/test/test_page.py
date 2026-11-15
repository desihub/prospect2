#
# See top-level LICENSE.rst file for Copyright information

"""
Tier 2: make_payload / make_html on synthetic Spectra.  Models are passed in and
alternate-fit models are off, so no redrock templates are needed.
"""

import json
import numpy as np
import pytest

from desispec.coaddition import coadd_cameras
from desispec.spectra import stack

from prospect2.page import make_payload, make_html, _fit_table
from prospect2.test.util import (synthetic_spectra, synthetic_zcat, synthetic_models, synthetic_fits,
                                 decode_q2, payload_from_html, TARGETID0)

NOALT = dict(alt_models=0)


@pytest.fixture(scope='module')
def sp():
    return synthetic_spectra(nspec=4)


def test_payload_basics(sp):
    zcat = synthetic_zcat(sp)
    p = make_payload(sp, zcat=zcat, models=synthetic_models(sp), provenance=dict(source='coadd-x.fits'), **NOALT)
    assert p['nspec'] == 4 and p['title'] == 'coadd-x.fits' and p['source'] == 'coadd-x.fits'
    assert p['model_source'] == 'input models'
    #- TARGETIDs > 2**53 go to JS as exact strings
    assert [t['targetid'] for t in p['targets']] == [str(TARGETID0 + i) for i in range(4)]
    assert all(t['z'] == 0.5 and t['spectype'] == 'GALAXY' and t['masked_cams'] == [] for t in p['targets'])
    assert [o[2] for o in p['overlaps']] == ['b/r overlap', 'r/z overlap']
    json.dumps(p)          #- serializable


def test_source_never_null(sp):
    #- the JS keys stored VI and names the VI export by `source`
    p = make_payload(sp, **NOALT)
    assert p['title'] == 'prospect2' and p['source'] == 'prospect2'
    p = make_payload(sp, title='My page', **NOALT)
    assert p['source'] == 'My page'
    p = make_payload(sp, title='My page', provenance=dict(source='coadd-x.fits'), **NOALT)
    assert p['title'] == 'My page' and p['source'] == 'coadd-x.fits'


def test_minimal_fibermap(sp):
    #- only TARGETID is required; other columns fall back to None / '' / defaults
    from astropy.table import Table
    sp2 = sp[:]
    sp2.fibermap = Table(dict(TARGETID=sp.fibermap['TARGETID']))
    p = make_payload(sp2, **NOALT)
    t = p['targets'][0]
    assert t['targetid'] == str(TARGETID0) and t['fiber'] is None and t['tileid'] is None
    assert t['ra'] is None and t['dec'] is None and t['ebv'] is None and t['exptime'] is None
    assert t['objtype'] == '' and t['morph'] == '' and t['mag']['r'] is None and t['fibermag_r'] is None
    assert t['numexp'] == 1 and t['fiberstatus'] == 0
    json.dumps(p)
    #- with a zcat and models too
    p = make_payload(sp2, zcat=synthetic_zcat(sp2), models=synthetic_models(sp2), **NOALT)
    assert p['targets'][0]['z'] == 0.5


def test_uncoadded_fibermap(sp):
    #- per-exposure spectra: EXPTIME / FIBERSTATUS instead of COADD_*, a single exposure
    sp2 = sp[:]
    fm = sp2.fibermap.copy()
    fm['EXPTIME'] = 900.0
    fm['FIBERSTATUS'] = 0
    fm['FIBERSTATUS'][1] = 2
    fm.remove_columns(['COADD_EXPTIME', 'COADD_NUMEXP', 'COADD_FIBERSTATUS'])
    sp2.fibermap = fm
    t = make_payload(sp2, **NOALT)['targets']
    assert t[0]['exptime'] == 900.0 and t[0]['numexp'] == 1 and t[0]['fiberstatus'] == 0
    assert t[1]['fiberstatus'] == 2


def test_targetid_required(sp):
    sp2 = sp[:]
    sp2.fibermap = sp.fibermap.copy()
    sp2.fibermap.remove_column('TARGETID')
    with pytest.raises(ValueError, match='TARGETID'):
        make_payload(sp2, **NOALT)


def test_payload_data_roundtrip(sp):
    p = make_payload(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp, 2.0), **NOALT)
    flux, noise, model, _ = decode_q2(p)
    merged = coadd_cameras(sp)
    f = merged.flux['brz']
    sigma = 1 / np.sqrt(merged.ivar['brz'])
    assert np.all(np.abs(flux - f) <= 0.051 * sigma)
    assert np.allclose(model, 2.0, atol=0.01)


def test_no_zcat(sp):
    p = make_payload(sp, **NOALT)
    assert p['model'] is None and p['redrock'] is None and 'z' not in p['targets'][0]


def test_spectra_attributes_and_overrides(sp):
    zcat = synthetic_zcat(sp)
    sp2 = sp[:]
    sp2.redshifts = zcat
    sp2.model = synthetic_models(sp, 3.0)
    p = make_payload(sp2, **NOALT)
    assert p['targets'][0]['z'] == 0.5
    assert np.allclose(decode_q2(p)[2], 3.0, atol=0.01)
    #- explicit arguments win
    zcat2 = synthetic_zcat(sp, z=1.25)
    p = make_payload(sp2, zcat=zcat2, models=synthetic_models(sp, 7.0), **NOALT)
    assert p['targets'][0]['z'] == 1.25
    assert np.allclose(decode_q2(p)[2], 7.0, atol=0.05)


def test_row_matching_errors(sp):
    zcat = synthetic_zcat(sp)
    with pytest.raises(ValueError, match='row-matched'):
        make_payload(sp, zcat=zcat[::-1], **NOALT)
    with pytest.raises(ValueError, match='row-matched'):
        make_payload(sp, zcat=zcat[:3], **NOALT)
    with pytest.raises(ValueError, match='shape'):
        make_payload(sp, zcat=zcat, models={b: m[:3] for b, m in synthetic_models(sp).items()}, **NOALT)
    with pytest.raises(ValueError, match='without zcat'):
        make_payload(sp, models=synthetic_models(sp), **NOALT)
    with pytest.raises(ValueError, match='entries'):
        make_payload(sp, zcat=zcat, models=synthetic_models(sp), fits=[None] * 3, **NOALT)


def test_list_input(sp):
    zcat = synthetic_zcat(sp)
    a, b = sp[:2], sp[2:]
    za, zb = zcat[:2], zcat[2:]
    models = synthetic_models(sp)
    ma = {k: v[:2] for k, v in models.items()}
    mb = {k: v[2:] for k, v in models.items()}
    one = make_html(sp, zcat=zcat, models=models, **NOALT)
    assert make_html([a, b], zcat=[za, zb], models=[ma, mb], **NOALT) == one
    #- lists must match in length; template versions must agree
    with pytest.raises(ValueError, match='list of 2'):
        make_html([a, b], zcat=[za], **NOALT)
    zb2 = zb.copy()
    zb2.meta['TEMVER00'] = '9.9'
    with pytest.raises(ValueError, match="TEMVER00"):
        make_html([a, b], zcat=[za, zb2], models=[ma, mb], **NOALT)
    #- the same check applies to spectra.redshifts
    a.redshifts, b.redshifts = za, zb2
    with pytest.raises(ValueError, match="TEMVER00"):
        make_html([a, b], **NOALT)


def test_masked_spectra():
    sp = synthetic_spectra(nspec=3)
    for b in sp.bands:                    #- spectrum 0: no data at all
        sp.ivar[b][0] = 0
    sp.ivar['z'][1] = 0                   #- spectrum 1: z camera masked
    p = make_payload(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp), **NOALT)
    flux = decode_q2(p)[0]
    t0, t1, t2 = p['targets']
    assert np.all(np.isnan(flux[0]))
    assert [m[2] for m in t0['masked_cams']] == ['b', 'r', 'z'] and t0['snr'] == dict(b=None, r=None, z=None)
    assert [m[2] for m in t1['masked_cams']] == ['z']
    assert t2['masked_cams'] == []
    #- spectrum 1 still has data in the r/z overlap (from r), none beyond r
    assert np.all(np.isnan(flux[1][p['nwave'] - 10:]))


def test_fit_table():
    d = synthetic_fits([0.5, 0.5016, 1.2, 0.5005], [100.0, 104.5, 110.0, 120.0])
    rows = _fit_table(d)
    assert [r[0] for r in rows] == [0.5, 0.5016, 1.2, 0.5005]
    assert [r[3] for r in rows] == [0.0, 4.5, 10.0, 20.0]                 # chi2 - best
    #- dup: index of a better fit within 1000 km/s (0.5016 is ~320 km/s from 0.5)
    assert [r[5] for r in rows] == [-1, 0, -1, 0]


def test_alternate_fits_display():
    sp = synthetic_spectra(nspec=3)
    sp.fibermap['OBJTYPE'][1] = 'SKY'
    zcat = synthetic_zcat(sp)
    fits = [synthetic_fits([0.5, 1.2], [100.0, 110.0]),
            synthetic_fits([0.5, 1.2], [100.0, 110.0]),
            None]
    p = make_payload(sp, zcat=zcat, models=synthetic_models(sp), fits=fits, **NOALT)
    t0, t1, t2 = p['targets']
    assert len(t0['fits']) == 2 and t0['fits'][1][:2] == [1.2, 'GALAXY']
    assert 'fits' not in t1          #- sky fibers get no alternate fits
    assert 'fits' not in t2
    assert 'amodel' not in p


def test_html(sp):
    html = make_html(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp), title='My page', **NOALT)
    assert '<meta charset="utf-8">' in html.lower()
    assert '<title>My page · prospect</title>' in html
    for placeholder in ('/*__UPLOT_CSS__*/', '/*__UPLOT_JS__*/', '"__DATA__"', '__TITLE__'):
        assert placeholder not in html
    assert payload_from_html(html)['title'] == 'My page'


def test_html_escaping():
    """Titles and payload strings can't break the page's HTML or end its <script>"""
    sp = synthetic_spectra(nspec=2)
    sp.fibermap['MORPHTYPE'] = ['</script><!--', 'REX']      # (replaces the 4-char column)
    clean = make_html(sp, alt_models=0, title='x')
    nasty = make_html(sp, alt_models=0, title='a < b & "c" </title>')
    assert '<title>a &lt; b &amp; &quot;c&quot; &lt;/title&gt; · prospect</title>' in nasty
    #- no extra script ends or comment openers beyond those in the template itself
    for tag in ('</script', '<!--'):
        assert nasty.count(tag) == clean.count(tag)
    p = payload_from_html(nasty)
    assert p['title'] == 'a < b & "c" </title>'
    assert p['targets'][0]['morph'] == '</script><!--'

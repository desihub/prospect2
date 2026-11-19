#
# See top-level LICENSE.rst file for Copyright information

"""
Tier 2: PageData (prospect2.pagedata): selection, stacking, HTML round trips,
edits, refusals, and the payload key-set guard.  Synthetic inputs only.
"""

import json
import numpy as np
import pytest

from prospect2 import PageData, stack, read_page, write_page, make_html
from prospect2.pagedata import FORMAT_VERSION, check_readable
from prospect2.test.util import (synthetic_spectra, synthetic_zcat, synthetic_models, synthetic_fits,
                                 synthetic_page, decode_q2, payload_from_html, TARGETID0)

NOALT = dict(alt_models=0)


def _from(sp, **kw):
    return PageData.from_spectra(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp),
                                 provenance=dict(source='coadd-x.fits'), **NOALT, **kw)


def test_basics():
    sp = synthetic_spectra(nspec=4)
    p = _from(sp)
    assert len(p) == p.num_spectra() == 4 and p.payload()['format_version'] == FORMAT_VERSION
    assert list(p.fibermap['TARGETID']) == [TARGETID0 + i for i in range(4)]
    assert p.fibermap['TARGETID'].dtype == np.int64 and list(p.fibermap['FIBER']) == [0, 1, 2, 3]
    assert list(p.redshifts['SPECTYPE']) == ['GALAXY'] * 4 and 'TARGETID' in p.redshifts.colnames
    assert p.flux.shape == (4, p.nwave) and p.model.shape == (4, p.nwave) and len(p.wave) == p.nwave
    assert np.allclose(p.model, 1.0, atol=0.01)


def test_select_equals_direct():
    #- selecting rows of a page == a page made directly from those spectra (incl. reordering)
    sp = synthetic_spectra(nspec=5)
    p = _from(sp)
    for sel in ([3, 1], np.array([False, True, True, False, True]), slice(1, 4)):
        assert p[sel].payload() == _from(sp[sel]).payload()


def test_stack_equals_direct():
    sp = synthetic_spectra(nspec=5)
    p = _from(sp)
    assert stack([p[:2], p[2:]]).payload() == p.payload()
    #- from separately generated pages too, and sources merge when identical
    s = stack([_from(sp[:3]), _from(sp[3:])])
    assert s.payload() == p.payload() and len(s.sources) == 1


def test_sources():
    sp = synthetic_spectra(nspec=4)
    a = PageData.from_spectra(sp[:2], provenance=dict(source='a.fits'), **NOALT)
    b = PageData.from_spectra(sp[2:], provenance=dict(source='b.fits'), **NOALT)
    s = stack([a, b])
    assert [x['source'] for x in s.sources] == ['a.fits', 'b.fits'] and list(s.extra['src']) == [0, 0, 1, 1]
    d = s.payload()
    assert d['title'] == 'prospect2' and d['vi_key'] == 'prospect2' and d['targets']['src'] == [0, 0, 1, 1]
    #- title: settable, from stack's argument; it keys the VI only without a single source
    assert stack([a, b], title='Two').payload()['title'] == 'Two' == stack([a, b], title='Two').title
    assert stack([a, b], title='Two').payload()['vi_key'] == 'Two'
    one = stack([a[:1], a[1:]], title='Subset')
    assert one.payload()['title'] == 'Subset' and one.payload()['vi_key'] == 'a.fits'
    #- selecting drops unused sources and renumbers
    t = s[[3, 2]]
    assert [x['source'] for x in t.sources] == ['b.fits'] and list(t.extra['src']) == [0, 0]
    assert t.payload()['vi_key'] == 'b.fits'


def test_alt_model_rows():
    p = synthetic_page(nspec=4)              #- alt rows for (0, 1), (2, 1), (2, 2)
    amodel = lambda q: np.cumsum(q.arrays['amodel']['i8'].astype(int), axis=1)   #- no escapes for these ramps
    q = p[[2, 1, 0]]
    assert q.arrays['amodel']['owner'].tolist() == [[0, 1], [0, 2], [2, 1]]
    assert np.array_equal(amodel(q), amodel(p)[[1, 2, 0]])
    assert np.array_equal(q.arrays['amodel']['scale'], p.arrays['amodel']['scale'][[1, 2, 0]])
    assert len(p[[1, 3]].arrays['amodel']['owner']) == 0
    #- stacking offsets owners; the result equals the original
    s = stack([p[:1], p[1:]])
    assert s.arrays['amodel']['owner'].tolist() == [[0, 1], [2, 1], [2, 2]]
    assert s.payload() == p.payload()


def test_html_roundtrip(tmp_path):
    p = synthetic_page(nspec=4)
    q = PageData.from_html(p.to_html())
    assert q.payload() == p.payload()
    assert np.array_equal(q.flux, p.flux, equal_nan=True) and q.title == 'coadd-x.fits'
    #- read_page / write_page
    write_page(str(tmp_path / 'sub' / 'page.html'), p[[0, 2]])
    r = read_page(str(tmp_path / 'sub' / 'page.html'))
    assert list(r.fibermap['FIBER']) == [0, 2] and r.payload()['arrays'] == p[[0, 2]].payload()['arrays']


def test_make_html_is_pagedata():
    sp = synthetic_spectra(nspec=3)
    zcat, models = synthetic_zcat(sp), synthetic_models(sp)
    html = make_html(sp, zcat=zcat, models=models, title='T', **NOALT)
    assert html == PageData.from_spectra(sp, zcat=zcat, models=models, **NOALT).to_html(title='T')


def test_column_roundtrip():
    #- dtypes survive: int64 > 2**53 exact, float32 exact, masked ints/strings, NaN floats
    sp = synthetic_spectra(nspec=3)
    p = _from(sp)
    p.fibermap['DESI_TARGET'] = np.array([2**62 + 3, 0, 4], dtype=np.int64)
    p.fibermap['FLUX_G'] = np.array([0.1, np.nan, 3e-5], dtype=np.float32)
    p.fibermap['OBJTYPE'] = np.ma.masked_array(['TGT', 'SKY', ''], mask=[False, False, True])
    q = PageData.from_payload(json.loads(json.dumps(p.payload())))
    assert q.fibermap['DESI_TARGET'].dtype == np.int64 and q.fibermap['DESI_TARGET'][0] == 2**62 + 3
    assert q.fibermap['FLUX_G'].dtype == np.float32
    assert np.array_equal(q.fibermap['FLUX_G'], p.fibermap['FLUX_G'], equal_nan=True)
    assert q.fibermap['OBJTYPE'].mask.tolist() == [False, False, True] and q.fibermap['OBJTYPE'][1] == 'SKY'
    for c in p.fibermap.colnames:
        assert q.fibermap[c].dtype.kind == p.fibermap[c].dtype.kind, c


def test_edits_are_written():
    #- derived values (target classes, ZWARN names) are recomputed from the tables
    from desitarget.targetmask import desi_mask
    p = _from(synthetic_spectra(nspec=2))
    p.fibermap['DESI_TARGET'] = [desi_mask.QSO, desi_mask.ELG]
    p.redshifts['ZWARN'][1] = 2           #- LITTLE_COVERAGE
    d = p.payload()
    assert d['targets']['classes'] == [['QSO'], ['ELG']] and d['targets']['zwarn_names'][1] == ['LITTLE_COVERAGE']
    p.fibermap['vi_flag'] = 1
    with pytest.raises(ValueError, match='UPPERCASE'):
        p.payload()


def test_stack_optional_parts():
    #- a page without top fits / alternate-fit models stacks with one that has them
    a = synthetic_page(nspec=3, owner=((0, 1), (1, 2)))
    sp = synthetic_spectra(nspec=2)
    sp.fibermap['FIBER'] = [10, 11]
    b = _from(sp)
    s = stack([a, b])
    assert len(s) == 5 and list(s.extra['zfits'][3:]) == [None, None]
    assert s.arrays['amodel']['owner'].tolist() == a.arrays['amodel']['owner'].tolist()
    assert s.alt_models == a.alt_models
    PageData.from_html(s.to_html())


def test_stack_refusals():
    sp = synthetic_spectra(nspec=4)
    p = _from(sp)
    with pytest.raises(ValueError, match='more than once'):
        stack([p, p])
    with pytest.raises(ValueError, match='no model'):     #- no zcat: no model, no redshifts
        stack([p[:2], PageData.from_spectra(sp[2:], **NOALT)])
    q = _from(sp[2:])
    q.wave0 += 1
    with pytest.raises(ValueError, match='wavelength grid'):
        stack([p[:2], q])
    q = _from(sp[2:])
    q.spgrp = 'healpix'
    p2 = p[:2]
    p2.spgrp = 'cumulative'
    with pytest.raises(ValueError, match='SPGRP'):
        stack([p2, q])
    with pytest.raises(ValueError, match='no pages'):
        stack([])


def test_version_checks():
    check_readable(FORMAT_VERSION)
    major, minor = map(int, FORMAT_VERSION.split('.'))
    assert check_readable(f'{major}.{minor + 1}') is True and check_readable(FORMAT_VERSION) is False
    for v, msg in ((None, 'no format_version'), (f'{major + 1}.0', 'upgrade'), ('x', 'malformed')):
        with pytest.raises(ValueError, match=msg):
            check_readable(v)
    d = _from(synthetic_spectra(nspec=2)).payload()
    with pytest.raises(ValueError, match='upgrade'):
        PageData.from_payload(dict(d, format_version=f'{major + 1}.0'))
    #- a stored column must be in the columns map
    with pytest.raises(ValueError, match='columns map'):
        PageData.from_payload(dict(d, targets=dict(d['targets'], MYSTERY=[0, 0])))
    with pytest.raises(ValueError, match='prospect2-data'):
        PageData.from_html('<html>const D = {}</html>')


def test_reordered_tables_refused():
    #- sorting a table in place would pair spectra with other targets' rows: refused when written
    p = _from(synthetic_spectra(nspec=3))
    p.redshifts['Z'] = [0.9, 0.1, 0.5]
    p.redshifts.sort('Z')
    with pytest.raises(ValueError, match='no longer match'):
        p.payload()
    with pytest.raises(ValueError, match='no longer match'):
        stack([p])
    q = _from(synthetic_spectra(nspec=3))
    q.redshifts['Z'] = [0.9, 0.1, 0.5]
    r = q[np.argsort(q.redshifts['Z'])]          #- the supported way
    assert list(r.redshifts['Z']) == [0.1, 0.5, 0.9] and r.payload()['targets']['Z'] == [0.1, 0.5, 0.9]


def test_other_row_space():
    #- a non-alt row space (e.g. future archetype models): owner column 0 = spectrum, other columns
    #- are carried as is and not checked against zfits; slicing and stacking use only column 0
    def with_arch(p, owner):
        a = dict(p.arrays['amodel'], rows='arch', owner=np.array(owner, dtype=np.int64))
        return PageData(p.fibermap, p.redshifts, p.extra, dict(p.arrays, arch=a), p.wave0, p.dwave, p.nwave,
                        p.sources)
    p = synthetic_page(nspec=4)
    q = with_arch(p, [[0, 7, 5], [2, 9, 5], [2, 8, 6]])      #- 7, 8, 9: not fit indices
    r = q[[3, 2, 0]]
    assert r.arrays['arch']['owner'].tolist() == [[1, 9, 5], [1, 8, 6], [2, 7, 5]]
    s = stack([q[:2], q[2:]])
    assert s.arrays['arch']['owner'].tolist() == q.arrays['arch']['owner'].tolist()
    assert np.array_equal(s.arrays['arch']['i8'], q.arrays['arch']['i8'])
    with pytest.raises(ValueError, match='spectrum index'):
        with_arch(p, [[4, 0, 0], [0, 0, 0], [0, 0, 0]])


def test_pr_review_fixes():
    #- stack leaves its inputs unchanged (owner offsets were applied in place)
    a, b = synthetic_page(nspec=4), synthetic_page(nspec=4, FIBER=[10, 11, 12, 13])
    before = b.arrays['amodel']['owner'].copy()
    s1 = stack([a, b])
    assert np.array_equal(b.arrays['amodel']['owner'], before)
    assert stack([a, b]).payload() == s1.payload()
    #- selecting from a page whose fibermap was sorted in place is refused (it would adopt the sorted order)
    p = synthetic_page(nspec=4)
    p.redshifts = None
    p.fibermap.sort('FIBER', reverse=True)
    with pytest.raises(ValueError, match='no longer match'):
        p[:]
    #- ROW VI keys hash TARGETID + grouping values: two observations of one target in two orders differ
    p = synthetic_page(nspec=2, owner=(), TARGETID=[5, 5], NIGHT=[20240101, 20240102])
    k1, k2 = p[[0, 0, 1]].payload()['vi_key'], p[[1, 0, 0]].payload()['vi_key']
    assert ':rows-' in k1 and k1 != k2


def test_fibermap_value_masked_fallback():
    #- after a stack with an outer join, a masked first-choice cell falls through to the next column
    from astropy.table import Table, MaskedColumn
    from prospect2.util import fibermap_value
    fm = Table([MaskedColumn([900.0, 0.0], name='COADD_EXPTIME', mask=[False, True]),
                MaskedColumn([0.0, 300.0], name='EXPTIME', mask=[True, False])])
    assert [fibermap_value(fm, i, 'COADD_EXPTIME', 'EXPTIME') for i in range(2)] == [900.0, 300.0]


def test_validation():
    p = synthetic_page(nspec=4)
    a = dict(p.arrays['amodel'], owner=np.array([[3, 1], [0, 1], [2, 1]]))     #- spectrum 3 has no fits
    with pytest.raises(ValueError, match='owner'):
        PageData(p.fibermap, p.redshifts, p.extra, dict(p.arrays, amodel=a), p.wave0, p.dwave, p.nwave, p.sources)
    extra = p.extra.copy()
    extra['src'] = 1
    with pytest.raises(ValueError, match='sources'):
        PageData(p.fibermap, p.redshifts, extra, p.arrays, p.wave0, p.dwave, p.nwave, p.sources)
    #- selections don't share mutable metadata with their parent
    q = p[:2]
    q.sources[0]['source'] = 'other.fits'
    q.alt_models['n'] = 9
    assert p.sources[0]['source'] == 'coadd-x.fits' and p.alt_models['n'] == 2


def test_alt_model_settings_in_stack():
    #- pages that show alternate fits must agree on alternate-fit models (None = none stored)
    a = synthetic_page(nspec=3, owner=((0, 1),))
    b = synthetic_page(nspec=3, owner=(), FIBER=[10, 11, 12])          #- zfits but no models stored
    assert b.alt_models is None
    with pytest.raises(ValueError, match='alternate-fit model settings'):
        stack([a, b])
    #- a page without any zfits is compatible with either
    c = _from(synthetic_spectra(nspec=2))
    c.fibermap['FIBER'] = [20, 21]
    assert stack([a, c]).alt_models == a.alt_models and stack([b, c]).alt_models is None
    #- from FITS: the settings are kept even when no fit qualified for a model
    sp = synthetic_spectra(nspec=2)
    zfits = [synthetic_fits([0.5, 1.2], [100.0, 300.0])] * 2              #- Δχ² 200 > 100: no models
    d = PageData.from_spectra(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp), zfits=zfits)
    assert 'amodel' not in d.arrays and d.alt_models == dict(n=1, max_dchi2=100.0)


def test_empty_list_selection():
    p = _from(synthetic_spectra(nspec=3))
    assert len(p[[]]) == 0 and len(p[[i for i in range(3) if i > 5]]) == 0


def test_string_width_and_int64_roundtrip():
    p = _from(synthetic_spectra(nspec=2))
    q = PageData.from_html(p[p.redshifts['SPECTYPE'] == 'GALAXY'][:1].to_html())
    assert q.redshifts['SPECTYPE'].dtype == p.redshifts['SPECTYPE'].dtype
    #- any int64 beyond 2**53 goes to JS as a string, exactly
    p.fibermap['MYID'] = np.array([2**53 + 1, 5], dtype=np.int64)
    d = p.payload()
    assert d['targets']['MYID'] == [str(2**53 + 1), '5']
    assert PageData.from_payload(d).fibermap['MYID'][0] == 2**53 + 1


def test_newer_derived_keys_ignored():
    #- derived keys are added without a version bump, so a reader must ignore ones it doesn't know
    d = _from(synthetic_spectra(nspec=2)).payload()
    d2 = dict(d, new_summary={'x': 1}, targets=dict(d['targets'], mag_g=[20.1, 20.2], odd=[1]))
    assert PageData.from_payload(d2).payload() == d


def test_newer_minor_read_with_drops(caplog):
    #- a newer minor of the same major is read; keys this code doesn't know are dropped with a warning
    #- naming them, and the rewritten page is the current format without them
    import logging
    major, minor = map(int, FORMAT_VERSION.split('.'))
    d = _from(synthetic_spectra(nspec=2)).payload()
    targets = dict(d['targets'], NEWCOL=[1, 2], newstored=['a', 'b'])
    arrays = dict(d['arrays'], sky=dict(d['arrays']['flux'], enc='q2-sky'),
                  flux=dict(d['arrays']['flux'], width=d['nwave']))
    d2 = dict(d, format_version=f'{major}.{minor + 1}', vi_settings={'x': 1},
              columns=dict(d['columns'], photometry={'NEWCOL': 'i8'}), targets=targets, arrays=arrays,
              sources=[dict(s, catalog='zall.fits') for s in d['sources']])
    caplog.set_level(logging.WARNING)
    p = PageData.from_payload(d2)
    for k in ('newer than', 'vi_settings', 'columns table photometry', 'sources field catalog', 'array sky',
              'array flux field width', 'per-target newstored'):
        assert k in caplog.text, k
    assert p.payload() == d
    #- the input dict is not modified
    assert 'photometry' in d2['columns'] and 'sky' in d2['arrays'] and 'catalog' in d2['sources'][0]


def test_no_nan_in_payload():
    #- JSON.parse rejects NaN/Infinity, which would leave a blank page
    sp = synthetic_spectra(nspec=3)
    zfits = [synthetic_fits([0.5, 1.2], [100.0, np.inf])] * 3
    p = PageData.from_spectra(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp), zfits=zfits,
                              alt_models=0)
    p.alt_models = dict(n=2, max_dchi2=None)
    p.extra['snr'][0] = dict(b=np.nan, r=np.inf, z=1.0)
    html = p.to_html()
    d = payload_from_html(html)
    assert d['targets']['zfits'][0][1][3] is None and d['targets']['snr'][0] == dict(b=None, r=None, z=1.0)
    start = html.index('id="prospect2-data">')
    raw = html[start:html.index('</script>', start)]
    for bad in ('NaN', 'Infinity'):
        assert bad not in raw


def test_infinite_alt_max_dchi2():
    sp = synthetic_spectra(nspec=2)
    p = PageData.from_spectra(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp), alt_max_dchi2=np.inf,
                              zfits=[None, None])
    assert p.alt_models is None        #- no alternate fits, so no settings; with them: max_dchi2 None
    json.loads(json.dumps(p.payload(), allow_nan=False))


def test_empty_and_repeated_selections(tmp_path):
    p = synthetic_page(nspec=4)
    empty = p[np.zeros(4, dtype=bool)]
    assert len(empty) == 0 and len(empty.arrays['amodel']['owner']) == 0
    with pytest.raises(ValueError, match='no spectra'):
        empty.to_html()
    with pytest.raises(ValueError, match='no spectra'):
        write_page(str(tmp_path / 'newdir' / 'e.html'), empty)
    assert not (tmp_path / 'newdir').exists()
    #- an empty selection (e.g. a petal with no matches) stacks fine
    s = stack([empty, p[[1, 2]], empty])
    assert s.payload() == p[[1, 2]].payload()
    #- repeated indices: rows and alternate-fit rows are repeated; VI is then keyed by row, with
    #- a vi_key specific to this list of spectra
    r = p[[2, 2]]
    assert r.arrays['amodel']['owner'].tolist() == [[0, 1], [0, 2], [1, 1], [1, 2]]
    d = r.payload()
    assert d['keycols'] == ['ROW'] and d['vi_key'].startswith('coadd-x.fits:rows-')
    assert d['vi_key'] != p[[1, 1]].payload()['vi_key'] and p.payload()['vi_key'] == 'coadd-x.fits'


def test_stack_duplicates():
    from prospect2.test.test_grouping import tile_spectra
    a = PageData.from_spectra(tile_spectra(1496), **NOALT)
    with pytest.raises(ValueError, match='more than once'):
        stack([a, a[:1]])
    #- pages that already had duplicates (VI keyed by row) still stack, but not with themselves
    aa = PageData.from_spectra([tile_spectra(1496), tile_spectra(1496)], **NOALT)
    b = PageData.from_spectra(tile_spectra(2000), **NOALT)
    assert len(stack([aa, b])) == len(aa) + len(b)
    with pytest.raises(ValueError, match='more than once'):
        stack([aa, aa])
    #- the same tile with and without a known NIGHT are different spectra here (no false alarm)
    sp = tile_spectra(1496)
    sp.meta = dict(SPGRP='cumulative', TILEID=1496)
    c = PageData.from_spectra(sp, **NOALT)
    assert 'NIGHT' not in c.fibermap.colnames
    assert len(stack([a, c])) == len(a) + len(c)


#- Stored key sets of the current payload format.  If a payload change fails this test, bump
#- pagedata.FORMAT_VERSION (major: a key removed/renamed/changed in meaning or encoding; minor: keys
#- added, which readers must tolerate being missing), record the new sets here, and update
#- doc/implementation.md "Payload contents".
KEYSETS = {
    '1.0': dict(
        top={'format_version', 'generator', 'title', 'vi_key', 'sources', 'spgrp', 'keycols', 'wave0', 'dwave',
             'nwave', 'nspec', 'flux_unit', 'columns', 'targets', 'lines', 'overlaps', 'alt_models', 'arrays'},
        source={'source', 'redrock', 'model_source', 'rrdetails'},
        fibermap={'TARGETID', 'SURVEY', 'PROGRAM', 'TILEID', 'FIBER', 'NIGHT', 'EXPID', 'TARGET_RA', 'TARGET_DEC', 'OBJTYPE', 'MORPHTYPE',
                  'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_R', 'EBV', 'COADD_EXPTIME',
                  'EXPTIME', 'COADD_NUMEXP', 'COADD_FIBERSTATUS', 'FIBERSTATUS',
                  'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',
                  'SV1_DESI_TARGET', 'SV3_BGS_TARGET', 'CMX_TARGET'},
        redshifts={'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'DELTACHI2', 'CHI2', 'NPIXELS'},
        lowercase={'snr', 'masked_cams', 'zfits', 'src', 'classes', 'primary', 'zwarn_names', 'fiberstatus',
                   'fiberstatus_names', 'exptime', 'numexp'},
        arrays={'noise': {'rows', 'pack', 'enc', 'data'}, 'flux': {'rows', 'pack', 'enc', 'data'},
                'model': {'rows', 'pack', 'enc', 'data', 'scale'},
                'amodel': {'rows', 'pack', 'enc', 'data', 'scale', 'owner'}},
        zfits_row=6,
    ),
}


def test_payload_keyset_guard():
    assert FORMAT_VERSION in KEYSETS, f'record the key sets of format {FORMAT_VERSION} in KEYSETS'
    k = KEYSETS[FORMAT_VERSION]
    #- every optional part present, and the optional grouping and target-bit columns
    d = synthetic_page(nspec=4, NIGHT=20240214, EXPID=[1, 2, 3, 4], PROGRAM='dark', SV1_DESI_TARGET=0,
                       SV3_BGS_TARGET=0, CMX_TARGET=0).payload()
    assert set(d) == k['top']
    assert all(set(s) == k['source'] for s in d['sources'])
    #- synthetic_spectra's fibermap is desispec's full blank one (plus the columns above), so this
    #- is the stored subset
    assert set(d['columns']['fibermap']) == k['fibermap']
    assert set(d['columns']['redshifts']) == k['redshifts']
    assert {c for c in d['targets'] if c != c.upper()} == k['lowercase']
    assert {n: set(a) for n, a in d['arrays'].items()} == k['arrays']
    assert {len(f) for z in d['targets']['zfits'] if z for f in z} == {k['zfits_row']}

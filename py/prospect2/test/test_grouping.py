#
# See top-level LICENSE.rst file for Copyright information

"""
Tier 2: spectral grouping (SPGRP) detection and the unique key that keys VI.
"""

import numpy as np
import pytest
from astropy.table import Table

from prospect2.grouping import grouping_columns, key_columns
from prospect2.page import make_payload, make_html
from prospect2.util import select_spectra
from prospect2.test.util import synthetic_spectra, payload_from_html, target_rows, TARGETID0

NOALT = dict(alt_models=0)


def tile_spectra(tileid, night=20240214, nspec=3):
    """synthetic cumulative tile coadd: header SPGRP/NIGHT, fibermap header FAPRGRM"""
    sp = synthetic_spectra(nspec=nspec)
    sp.fibermap['TILEID'] = tileid
    sp.fibermap.meta.update(SURVEY='main', PROGRAM='DARK1B', FAPRGRM='dark1b')
    sp.meta = dict(SPGRP='cumulative', SPGRPVAL=night, NIGHT=night, TILEID=tileid)
    return sp


def healpix_spectra(survey, program, nspec=3):
    """synthetic uniqpix coadd: no FIBER/TILEID columns, SURVEY/PROGRAM only in the header"""
    sp = synthetic_spectra(nspec=nspec)
    sp.fibermap.remove_columns(['FIBER', 'TILEID'])
    sp.fibermap.meta.update(SURVEY=survey, PROGRAM='DARK1B')
    sp.meta = dict(SPGRP='uniqpix', SPGRPVAL=9901, SURVEY=survey, PROGRAM=program)
    return sp


def test_tile_grouping():
    sp = tile_spectra(1496)
    g = grouping_columns(sp)
    assert g['SPGRP'][0] == 'cumulative' and g['NIGHT'] == [20240214] * 3
    assert g['SURVEY'][0] == 'main' and g['PROGRAM'][0] == 'dark1b'     #- the tile's FAPRGRM
    assert key_columns(g, sp.fibermap['TARGETID']) == ('cumulative', ['TILEID', 'FIBER', 'NIGHT'])
    #- NIGHT from SPGRPVAL when there is no NIGHT keyword
    del sp.meta['NIGHT']
    assert grouping_columns(sp)['NIGHT'][0] == 20240214


def test_healpix_grouping():
    sp = healpix_spectra('main', 'dark')
    g = grouping_columns(sp)
    assert g['PROGRAM'][0] == 'dark' and g['FIBER'][0] is None    #- primary header, not fibermap's
    assert key_columns(g, sp.fibermap['TARGETID']) == ('uniqpix', ['SURVEY', 'PROGRAM', 'TARGETID'])


def test_uncoadded_grouping():
    #- per-exposure rows: the same fiber repeats, NIGHT/EXPID columns make it unique
    sp = tile_spectra(1496, nspec=4)
    sp.fibermap['FIBER'] = [10, 10, 11, 11]
    sp.fibermap['NIGHT'] = 20240214
    sp.fibermap['EXPID'] = [100, 101, 100, 101]
    g = grouping_columns(sp)
    assert key_columns(g, sp.fibermap['TARGETID'])[1] == ['TILEID', 'FIBER', 'NIGHT', 'EXPID']


def test_unknown_grouping():
    #- no SPGRP: TILEID+FIBER if known, else TARGETID; ROW (warned) if neither is unique
    sp = synthetic_spectra(nspec=3)
    assert key_columns(grouping_columns(sp), sp.fibermap['TARGETID'])[1] == ['TILEID', 'FIBER']
    sp.fibermap = Table(dict(TARGETID=sp.fibermap['TARGETID']))
    assert key_columns(grouping_columns(sp), sp.fibermap['TARGETID']) == (None, ['TARGETID'])
    assert key_columns(grouping_columns(sp), [1, 1, 2])[1] == ['ROW']


def test_combined_tiles():
    #- the same targets on two tiles: one key each, read from each input's header before stacking
    a, b = tile_spectra(1496), tile_spectra(2000, night=20250101)
    p = payload_from_html(make_html([a, b], **NOALT))
    assert p['spgrp'] == 'cumulative' and p['keycols'] == ['TILEID', 'FIBER', 'NIGHT']
    t = target_rows(p)
    assert [(x['TILEID'], x['NIGHT']) for x in t] == [(1496, 20240214)] * 3 + [(2000, 20250101)] * 3
    assert t[0]['TARGETID'] == t[3]['TARGETID'] == str(TARGETID0)
    assert t[0]['SURVEY'] == 'main' and t[0]['PROGRAM'] == 'dark1b'


def test_combined_healpix():
    #- the same TARGETID in two programs is two spectra
    a, b = healpix_spectra('main', 'dark'), healpix_spectra('main', 'bright')
    p = payload_from_html(make_html([a, b], **NOALT))
    assert p['keycols'] == ['SURVEY', 'PROGRAM', 'TARGETID']
    assert p['targets']['PROGRAM'] == ['dark'] * 3 + ['bright'] * 3
    assert 'FIBER' not in p['targets'] and 'TILEID' not in p['targets']


def test_duplicates_fall_back_to_row():
    a = tile_spectra(1496)
    p = payload_from_html(make_html([a, a], **NOALT))
    assert p['keycols'] == ['ROW']


def test_mixed_groupings():
    with pytest.raises(ValueError, match='mix spectral groupings'):
        make_html([tile_spectra(1496), healpix_spectra('main', 'dark')], **NOALT)


def test_prestacked_keeps_first_header():
    #- a pre-stacked Spectra has at most the first input's headers (documented in make_html):
    #- columns like TILEID survive, header-only values like NIGHT don't
    from desispec.spectra import stack
    sp = stack([tile_spectra(1496), tile_spectra(2000, night=20250101)])
    t = target_rows(make_payload(sp, **NOALT))
    assert t[3]['TILEID'] == 2000 and t[3].get('NIGHT') != 20250101


def test_select_fibers_needs_fiber_column():
    fm = healpix_spectra('main', 'dark').fibermap
    with pytest.raises(ValueError, match='no FIBER column'):
        select_spectra(fm, fibers='1-3')
    assert select_spectra(fm, targetids=str(TARGETID0)).sum() == 1

#
# See top-level LICENSE.rst file for Copyright information

"""Tier 1: fiber lists, spectrum selection, output naming"""

import os
import numpy as np
import pytest
from astropy.table import Table

from prospect2.util import parse_fiberlist, select_spectra, default_outfile, resolve_outfile


def test_parse_fiberlist():
    assert list(parse_fiberlist('10,20')) == [10, 20]
    assert list(parse_fiberlist('30-33')) == [30, 31, 32, 33]          # inclusive
    assert list(parse_fiberlist('50:53')) == [50, 51, 52]              # python-style
    assert list(parse_fiberlist(' 1, 3-4 ,6:8,')) == [1, 3, 4, 6, 7]


def test_select_spectra():
    fm = Table(dict(FIBER=[1, 2, 3, 4], TARGETID=[100, 200, 300, 400]))
    assert select_spectra(fm).all()
    assert list(select_spectra(fm, fibers='2-3')) == [False, True, True, False]
    assert list(select_spectra(fm, targetids='400')) == [False, False, False, True]
    #- -f and -t select the union
    assert list(select_spectra(fm, fibers='1', targetids='300,400')) == [True, False, True, True]
    assert not select_spectra(fm, fibers='99').any()


@pytest.mark.parametrize('infile, expected', [
    ('coadd-4-1496-thru20240214.fits', 'prospect-4-1496-thru20240214.html'),
    ('/prod/tiles/coadd-main-dark-10032.fits.gz', 'prospect-main-dark-10032.html'),
    ('spectra-x.fits', 'prospect-x.html'),
    ('myspectra.fits', None),
    ('coadd-.fits', None),
    ('coadd-x.h5', None),
])
def test_default_outfile(infile, expected):
    assert default_outfile(infile) == expected


def test_resolve_outfile(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    infile = '/prod/coadd-4-1.fits'
    assert resolve_outfile(None, infile, subset=False) == 'prospect-4-1.html'
    assert resolve_outfile('out/', infile, subset=False) == 'out/prospect-4-1.html'
    os.mkdir('existing')
    assert resolve_outfile('existing', infile, subset=False) == os.path.join('existing', 'prospect-4-1.html')
    #- explicit file names are used as given, with or without -f/-t
    assert resolve_outfile('x.html', infile, subset=True) == 'x.html'
    assert resolve_outfile('foo', infile, subset=False) == 'foo'       # (with a warning)


@pytest.mark.parametrize('outarg', [None, 'out/'])
def test_resolve_outfile_errors(outarg):
    #- a subset can't use the default name, nor can an input with another name
    with pytest.raises(ValueError, match='subset'):
        resolve_outfile(outarg, '/prod/coadd-4-1.fits', subset=True)
    with pytest.raises(ValueError, match="can't derive"):
        resolve_outfile(outarg, '/prod/myspectra.fits', subset=False)

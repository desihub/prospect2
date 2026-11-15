#
# See top-level LICENSE.rst file for Copyright information

"""Tier 2: plot_prospect2 command-line parsing and output naming"""

import os
import pytest

from prospect2.scripts.plot_prospect2 import parse


def test_parse_default_outfile(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = parse(['-i', '/prod/coadd-4-1496-thru20240214.fits'])
    assert args.outfile == 'prospect-4-1496-thru20240214.html'
    args = parse(['-i', '/prod/coadd-4-1496-thru20240214.fits', '-o', 'out/'])
    assert args.outfile == 'out/prospect-4-1496-thru20240214.html'
    args = parse(['-i', '/prod/coadd-x.fits', '-f', '1-3', '-o', 'page.html'])
    assert args.outfile == 'page.html' and args.fibers == '1-3'


@pytest.mark.parametrize('argv, message', [
    (['-i', 'coadd-x.fits', '-f', '1'], 'subset'),
    (['-i', 'coadd-x.fits', '-t', '123', '-o', 'out/'], 'subset'),
    (['-i', 'myspectra.fits'], "can't derive"),
])
def test_parse_errors(argv, message, capsys):
    with pytest.raises(SystemExit) as err:
        parse(argv)
    assert err.value.code == 2
    assert message in capsys.readouterr().err

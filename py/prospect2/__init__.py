"""
prospect2: offline HTML spectral viewer pages for DESI spectra.

Main entry points: prospect2.make_html(spectra, ...) -> HTML string (see
prospect2.page), and PageData / read_page / stack / write_page to load,
select, stack and write pages (see prospect2.pagedata).
"""

from ._version import __version__
from .page import make_html
from .pagedata import PageData, stack
from .io import read_page, write_page, write_html

"""
prospect2: offline HTML spectral viewer pages for DESI spectra.

Main entry point: prospect2.make_html(spectra, ...) -> HTML string;
see prospect2.page.
"""

__version__ = '0.1.dev'

from .page import make_html, write_html

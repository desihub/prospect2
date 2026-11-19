#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.page
==============

Build a prospect2 viewer page from in-memory objects.  The main entry point is
make_html(), which takes desispec Spectra (and optionally redrock outputs) and
returns the HTML of one self-contained page.  It is
PageData.from_spectra(...).to_html(); see prospect2.pagedata to select, stack
and rewrite pages.
"""

import json, html

from .io import DATADIR, write_html
from .pagedata import PageData, FORMAT_VERSION


def make_payload(spectra, zcat=None, models=None, zfits=None, title=None, provenance=None,
                 alt_models=2, alt_max_dchi2=100, model_precision=100, alt_model_precision=30):
    """
    The page's data payload (a JSON-serializable dict); see make_html for the
    arguments.  The keys are listed in doc/implementation.md "Payload
    contents"; keep that in sync.
    """
    return PageData.from_spectra(spectra, zcat=zcat, models=models, zfits=zfits, provenance=provenance,
                                 alt_models=alt_models, alt_max_dchi2=alt_max_dchi2,
                                 model_precision=model_precision,
                                 alt_model_precision=alt_model_precision).payload(title=title)


def render_html(payload):
    """
    Fill the viewer template with uPlot and `payload`; returns the page HTML.
    The payload goes in a <script type="application/json" id="prospect2-data">
    element.  The title is HTML-escaped, and every "<" in the payload JSON is
    written as \\u003c (same string in JSON), so no string can end the <script>
    element ("</script>") or open a comment ("<!--").  The data goes in last, so
    text in it is never mistaken for a placeholder.
    """
    template = (DATADIR / 'prospect2_template.html').read_text(encoding='utf-8')
    uplot_js = (DATADIR / 'vendor' / 'uPlot.iife.min.js').read_text(encoding='utf-8')
    uplot_css = (DATADIR / 'vendor' / 'uPlot.min.css').read_text(encoding='utf-8')
    #- allow_nan=False: JSON.parse rejects NaN/Infinity, which would leave a blank page
    data = json.dumps(payload, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')
    return (template.replace('__TITLE__', html.escape(payload['title']))
                    .replace('__FORMAT_VERSION__', FORMAT_VERSION)
                    .replace('/*__UPLOT_CSS__*/', uplot_css)
                    .replace('/*__UPLOT_JS__*/', uplot_js)
                    .replace('"__DATA__"', data))


def make_html(spectra, zcat=None, models=None, zfits=None, title=None, provenance=None,
              alt_models=2, alt_max_dchi2=100, model_precision=100, alt_model_precision=30):
    """
    HTML of a self-contained prospect2 viewer page.

    Args:
        spectra: desispec Spectra with per-camera data (not merged across
            cameras), e.g. from read_spectra(); or a list of them, which are
            stacked after checking that their redrock template versions agree.
            When passing one Spectra built by stacking several inputs yourself,
            you are responsible for that consistency, and only the first
            input's headers (SPGRP, NIGHT, SURVEY, PROGRAM) survive the
            stack, so pass the list instead.  Spectra are keyed for VI by
            TILEID/FIBER/NIGHT(/EXPID) for tile-based groupings (SPGRP
            cumulative, pernight, perexp) and by SURVEY/PROGRAM/TARGETID for
            healpix/uniqpix; mixing groupings raises ValueError.

    Options (each row-matched to spectra; for a list of spectra, a list of the same length):
        zcat: redrock REDSHIFTS table; default spectra.redshifts
            (read_spectra(..., return_redshifts=True)).  Its meta should hold the
            TEMNAMnn/TEMVERnn keywords so computed models use the fit's templates.
        models: per-camera best-fit models {band: [nspec, nwave_band]}, e.g. from
            prospect2.io.read_rrmodel; default spectra.model
            (read_spectra(..., return_models=True)); if neither, computed from
            the redrock templates + COEFF.
        zfits: the top redrock fits per spectrum, from prospect2.io.read_rrdetails.
        title: page title; default provenance['source'].
        provenance: file names shown on the page, keys source, redrock,
            model_source, rrdetails.
        alt_models: store models for alternate fits 2..alt_models+1 (0 = none).
        alt_max_dchi2: only for alternate fits with chi2 - best chi2 below this.
        model_precision, alt_model_precision: store models in steps of
            (median data noise)/N.

    Explicit zcat/models override spectra.redshifts/spectra.model.  Raises
    ValueError if the inputs are not row-matched or consistent.
    """
    page = PageData.from_spectra(spectra, zcat=zcat, models=models, zfits=zfits, provenance=provenance,
                                 alt_models=alt_models, alt_max_dchi2=alt_max_dchi2,
                                 model_precision=model_precision, alt_model_precision=alt_model_precision)
    return page.to_html(title=title)

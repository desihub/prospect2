# prospect2 prototype

Either install it (`pip install .`, or `pip install .[all]` to also get desitarget for target-bit
names), or use it in place with `bin/` on `$PATH` and `py/` on `$PYTHONPATH`:

```
plot_prospect2 -i coadd-4-1496-thru20240214.fits                 # -> ./prospect-4-1496-thru20240214.html
plot_prospect2 -i coadd-4-1496-thru20240214.fits -o outdir/       # -> outdir/prospect-4-1496-thru20240214.html
plot_prospect2 -i coadd-4-1496-thru20240214.fits -f 2000:2020 -o example-20.html
```

Matching `redrock-*.fits`, `rrmodel-*.fits` and `rrdetails-*.h5` files are found next to the coadd.
Use `-f`/`-t` to select a subset (see below).

From Python, e.g. to combine several coadds into one page:

```python
import prospect2
from desispec.io import read_spectra

spectra = [read_spectra(f, return_models=True, return_redshifts=True) for f in coadd_files]
spectra = [sp[sp.fibermap['OBJTYPE'] == 'TGT'] for sp in spectra]    # any selection
html = prospect2.make_html(spectra)                              # or a single Spectra
prospect2.write_html(html, 'mypage.html')
```

`make_html` takes a `Spectra` or a list of them (stacked after checking that the redrock template
versions agree), plus optional `zcat`, `models` and `zfits` (top redrock fits from `read_rrdetails`),
each row-matched to the spectra. `zcat` and `models` default to `spectra.redshifts` and
`spectra.model`; explicit arguments override them. Without models, they are computed from the
redrock templates. Only the flux and a fibermap with TARGETID are required: other fibermap
columns (coordinates, photometry, `COADD_*`, …) are shown as "—" when missing, and per-exposure
`EXPTIME`/`FIBERSTATUS` are used for uncoadded spectra. `provenance=dict(source=, redrock=,
model_source=, rrdetails=)` (all optional) names the input files; they are shown on hover over the
page name and in the `?` help.

### Filtering and stacking existing pages

Pages can be read back and filtered, combined and rewritten without the FITS files, in the same
way as desispec `Spectra` (`read_spectra` → select → `stack` → `write_spectra`). Reading a
500-spectrum page takes ~0.1 s, vs several seconds to generate it.

```python
import glob
from prospect2 import read_page, write_page, stack
from desitarget.targetmask import desi_mask

pages = [read_page(f) for f in sorted(glob.glob('prospect-*-1496-thru20240214.html'))]
qso = stack([p[(p.fibermap['DESI_TARGET'] & desi_mask.QSO) != 0] for p in pages])
write_page('qso-1496.html', qso)
hiz = qso[qso.redshifts['Z'] > 2]              # select like Spectra: mask, indices or slice
hiz.title = 'Tile 1496 QSOs, z > 2'            # label at the top left of the page
write_page('qso-1496-hiz.html', hiz)
```

A `PageData` has `fibermap` and `redshifts` (astropy Tables with DESI column names: TARGETID,
TILEID, FIBER, TARGET_RA, FLUX_G, DESI_TARGET, …, Z, SPECTYPE, ZWARN, …), `extra` (prospect2's
per-target values), `sources` (input files), and decoded `wave`, `flux`, `noise`, `model` arrays.
Selecting and stacking are lossless: a stacked page has exactly the data a page made from the same
spectra would. Edits to the tables are written; reorder only with `p[order]`, never by sorting a table in
place (that is caught only when it changes the TARGETID order). `stack` refuses pages with different wavelength
grids or groupings, model vs no model, or the same spectrum twice. The page title (top left) is
`stack(pages, title=...)`, `p.title = ...`, or `write_page(..., title=...)`; by default a selection
keeps its page's title, and a stack is named after its single source file, else the output file. Pages hold a data format
version (`prospect2.pagedata.FORMAT_VERSION`); pages from before format 1.0 must be regenerated. A page
from a newer prospect2 of the same major version is read, but data this version doesn't know are
dropped with a warning.

### Combining files

The grouping (header `SPGRP`) decides what identifies a spectrum, which keys the
VI and the VI export columns. Pass a *list* of Spectra, not one you stacked yourself: stacking keeps
only the first file's headers (SPGRP, NIGHT, SURVEY, PROGRAM).

| SPGRP                         | unique key                                  | list label                                |
|-------------------------------|---------------------------------------------|-------------------------------------------|
| cumulative, pernight, perexp  | TILEID, FIBER, NIGHT (+ EXPID if uncoadded) | fiber, or `tile:fiber` with several tiles |
| healpix, uniqpix              | SURVEY, PROGRAM, TARGETID                   | `#row`                                    |
| none (e.g. hand-made Spectra) | TILEID, FIBER if known, else TARGETID       | as above                                  |

Mixing groupings is an error. If the key is not unique (e.g. the same file twice), VI falls back to
the row on the page, with a warning. The same TARGETID elsewhere on the page (another tile, night or
program) is linked under the title ("also tile 2000 night 20250101").

Output name: `coadd-X.fits[.gz]` (or `spectra-X.fits[.gz]`) → `prospect-X.html`, written to the
current directory, or into `-o DIR` if that is an existing directory or ends in `/`. This default is
reserved for pages with every spectrum in the input, so with `-f`/`-t`, or an input named
differently, `-o FILE.html` is required. Any other `-o` is a file name (warning if it doesn't end in
`.html`/`.htm`). Missing directories are created; existing files are overwritten (and the log says
so).

Page sizes for petal 4, tile 1496:

| spectra | prospect (Bokeh)                 | prospect2           |
|---------|----------------------------------|---------------------|
| 20      | 4.95 MB (+ Bokeh JS from CDN)    | 0.51 MB, offline    |
| 500     | 107 MB (+ Bokeh JS from CDN)     | 8.87 MB, offline    |

Deep links: `#i=7&s=3` (target index, smoothing), `&view=grid`, `&theme=dark`, `&q=<filter>`.

Fiber/TARGETID lists, same syntax in `-f` (not for healpix coadds, which have no FIBER) and in the page's filter box: `10,20,30-40,50:60` =
fibers 10, 20, 30–40 (inclusive), 50–59 (python-style). In the filter box, integers ≥ 5000 or
negative are TARGETIDs (exact match, or partial match for a single number). `-f` and `-t` together
select the union.

Filter box: whitespace-separated terms, all of which must match (`2030-2040 ELG 0.79<z<0.81 !QSO`).
A term is a comparison on `z`, `dchi2`, `fiber`, `tile`, `night`, `snr` or `zwarn` (`z>1.2`,
`1.2<=z`, `0.79<z<0.81`, `tile=1496`), a fiber/TARGETID list (no spaces inside it), or a
case-insensitive text match on the target's labels (`ELG`, `GALAXY`, `nodata`, `unassigned`, survey
and program like `sv3`, `dark`). Prefix `!` negates a term. Terms must not contain spaces
(`z > 1.2` doesn't work). Terms that can't be parsed are shown in red under the box and ignored.

Dependencies: numpy, astropy, desispec, redrock (+ templates via `$RR_TEMPLATE_DIR`), and
optionally desitarget for target-bit names. No dependency on prospect: the line lists
(`py/prospect2/data/{emission,absorption}_lines.txt`) and `load_redrock_templates` were copied from it.
License: BSD 3-clause (`LICENSE.rst`), same as prospect.

## Tests

```
pytest                 # tiers 1-3 (about 5 s)
pytest -m browser      # tier 4: page in headless Chrome
```

Plain `pytest` never runs tier 4: it is opt-in (`addopts` in `pyproject.toml`), because it needs
Chrome and is slower. Run it separately after changing the viewer template or the payload
encoding; it is the only test of the JS.

| tier | file(s)                                      | needs                                      |
|------|----------------------------------------------|--------------------------------------------|
| 1    | `test_util.py`, `test_encode.py`, `test_models.py` | nothing                              |
| 2    | `test_page.py`, `test_pagedata.py`, `test_grouping.py`, `test_script.py` | nothing (synthetic Spectra) |
| 3    | `test_realdata.py`                           | petal 4 of tile 1496, else skipped         |
| 4    | `test_browser.py` (opt-in)                   | Chrome (`$CHROME`, the macOS app, or `google-chrome`/`chromium` on `$PATH`) |

Tier 3 looks for `coadd-4-1496-thru20240214.fits` and its redrock files in `$PROSPECT2_TEST_DATA`,
then `scratch/data/`, then `$DESI_ROOT/spectro/redux/{loa,nevis}/tiles/cumulative/1496/20240214/`.
Checks that need rrdetails or redrock templates (`$RR_TEMPLATE_DIR`) skip without them.

## Architecture

* Python (`prospect2.PageData`, used by `make_html` and `bin/plot_prospect2`) builds a versioned
  JSON payload, embedded as `<script type="application/json" id="prospect2-data">`;
  the HTML template (`py/prospect2/data/prospect2_template.html`) is a small hand-written JS app
  plus inlined uPlot (50 KB, `py/prospect2/data/vendor/`). No Bokeh, no server.
* Cameras merged with `desispec.coaddition.coadd_cameras`. Per-camera redrock models come from the
  matching `rrmodel*.fits` (`redrock-X.fits` → `rrmodel-X.fits`, rows matched by TARGETID, wavelength
  grids checked, error on mismatch) if present, else are computed from templates + COEFF (`camera_models`; forced with
  `--no-rrmodel`). Either way they're merged with the same ivar weights as the flux
  (`merge_camera_models`). The two agree except within ~5 px of masked pixels, where the rrmodel
  sags toward 0 (e.g. fiber 2183 z, pix 809: rrmodel 22.1, computed 29.6, data 29.4), as if the
  masked pixels were zeroed before the resolution convolution.
  prospect's `create_model` hard-cuts at mid-overlap, giving model=0 where that camera is masked
  but the other camera has good data (e.g. fiber 2208 at 7537.6 Å).
* Alternate redrock fits: if `rrdetails-X.h5` matches `redrock-X.fits`, the top N fits per target
  (`--nfits`, default 5; 0 ignores the file) are stored as z, spectype, subtype, Δχ² vs best, ZWARN
  (~170 bytes/spectrum). The page lists them under "Redrock fits" in the right sidebar (click a row, or `<` / `>` to step);
  selecting one sets z, so line markers and key-line panels follow it, and panels are chosen for its
  spectype. Fits within 1000 km/s of a better fit are faded and marked `≈#k`. Sky fibers get no
  alternate fits.
* Alternate-fit models: fits 2..N+1 (`--alt-models N`, default 2; 0 = none) with χ² − best χ² <
  `--alt-max-dchi2` (default 100) get their own redrock model, computed from the rrdetails
  coefficients + templates and stored at `--alt-model-precision` (default 30, i.e. σ/30). Sky and
  unassigned fibers never do. For other alternates the model is hidden and the plot says why.
  Defaults add ~3.5% to the 500-spectrum page (325 models).
* **Data encoding ("q2")**: noise as 1%-step log, delta-coded; flux quantized at 0.1σ; model
  delta-coded in steps of (median data σ)/N. All packed as int8 + int32 escapes, then gzip+base64;
  `--model-precision N` (default 100; span(model)/(10N) for spectra with no unmasked data).
  The model is noiseless, so steps can show as staircases: N=100 is invisible in the full view but
  can reach ~3–6 device px in zoomed, heavily smoothed flat regions; N=1000 stays < 0.4 px
  (+1.2 MB per 500 spectra). Flux error ≤ 0.05σ, mask preserved.
  ~10x smaller than the Bokeh page for the same 20 spectra.
* Trial-z model = pipeline model shifted in log λ (ignores resolution difference, fine for display).
* VI stored in localStorage, keyed by the grouping's unique key (above), exported as CSV with
  TARGETID, SURVEY, PROGRAM, then the rest of the key (TILEID, FIBER, NIGHT[, EXPID] for tiles).
* Line markers (`v` key / toolbar "Full | Data | Model | Edge", hash `&ls=`), default **Data**:
  SDSS-style short ticks, emission above / absorption below the plotted spectrum. The tick starts just
  outside everything plotted at the line (flux, model, and the noise trace above only); the label is
  also pushed clear of everything plotted under its full width, so labels never sit on noise spikes.
  Same-name lines within 60 px share a bracketed label (doublets); any lines within 10 px merge
  ("D2/D1", "K/H"). Colliding labels stack outward with longer leaders; peaks beyond the plot
  edge get a short stub. "Model" anchors ticks at model ± 2.5σ (steadier but crosses noise when
  unsmoothed); "Edge" puts ticks in strips at the plot top/bottom.
* Line-marker colors (`V` key, hash `&lc=`): **ink** (default; emission ink, absorption dark gray),
  plum (emission plum, absorption ink), crimson, violet (the original). The data are blue/orange, so
  any green or violet marker collides with them for colour-blind readers (green vs orange model
  ΔE≈1 protan; violet vs blue flux ΔE≈5 deutan); neutral ink is ≥25 ΔE from both in every CVD sim.
  Ticks/labels start 8 px from the data (`&lg=` px to experiment) and ticks carry a surface halo.
* Fully masked spectra (every pixel ivar=0, e.g. unassigned/broken positioners with COADD_NUMEXP=0):
  "NO DATA" badge, hatched main plot with the reason (NUMEXP, FIBERSTATUS), no key-line panels,
  hatched grid thumbnail. Redrock's placeholders for these (NPIXELS=0, CHI2=9e99, DELTACHI2≈1.9e84)
  are shown as "no fit" / "—" (raw values kept in the Redrock table). Blank OBJTYPE + negative
  TARGETID is labeled UNASSIGNED, not SKY. Filter tokens `nodata`, `unassigned`.
* Entirely masked camera (e.g. BADAMPZ: flux=ivar=0 for the whole z camera, 4 fibers in petal 4):
  the camera's range is hatched and labeled "z camera masked (BADAMPZ)", and key-line panels in it
  say "z camera masked". Neither data nor model are drawn there (the coadd stores zeros, and the
  PCA model extrapolated over a whole camera is distractingly wrong); nothing masked enters the
  y-range. Grid thumbnails hatch the span too. Small masked regions are still plain gaps.

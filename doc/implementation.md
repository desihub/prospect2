# prospect2 implementation notes

Key design decisions and the reasons behind them. `README.md` covers usage; this file covers *why*.
Keep it current when a decision changes, and record the reason (often a user preference or a
specific fiber that exposed a problem).

## Architecture

**Python builds a data payload; a static hand-written JS app renders it.** The generator
(`bin/plot_prospect2` → `prospect2.make_html`) reads the coadd/redrock files and
writes one self-contained HTML file: `prospect2_template.html` with uPlot 1.6.31 (`vendor/`,
~50 KB) and the payload inlined.

**Package layout follows desispec** (user decision, 2026-11-14): `bin/` wrappers without `.py`
that call `prospect2.scripts.X.parse()` / `main(args)`; library code in `py/prospect2/`; files the
code reads in `py/prospect2/data/`, found with `importlib.resources`, so it all works with
`bin/` on `$PATH` and `py/` on `$PYTHONPATH`, without installation. The viewer template
(`data/prospect2_template.html`) is really the app's source code, not data, but it lives in
`data/` with the other package files so installed copies can find it.

- *Why not patch Bokeh (the original prospect)?* The Bokeh page was 4.7 MB for 20 spectra plus
  CDN-loaded BokehJS, and its interaction model is hard to customize.
- *Why not Plotly?* The library alone is ~3.5 MB.
- *Why not Panel/Dash?* They need a running server; pages must work offline as a file, by email,
  or on a static web server.
- *Why uPlot?* It is tiny and fast for 500 × 7781-pixel spectra, and exposes canvas draw hooks for
  custom overlays (line markers, hatching).

**No dependency on prospect.** The line lists (`py/prospect2/data/{emission,absorption}_lines.txt`) and
`load_redrock_templates` were copied from prospect, so the prototype stands alone. The
dependencies are numpy, astropy, desispec, redrock (+ templates via `$RR_TEMPLATE_DIR`), and
optionally desitarget. License: BSD 3-clause, same authors as prospect.

**Packaging** (`pyproject.toml` only, modeled on desiutil, the most modern desihub layout; user
decision 2026-11-14 to add it before sharing):
- `bin/` scripts are installed with `script-files`, not `console_scripts` entry points, so the
  same `bin/` wrappers work installed and in place (DESI convention).
- Version is `__version__` in `py/prospect2/__init__.py`, read statically by setuptools (no import,
  so building doesn't need desispec).
- Dependencies: numpy, astropy, h5py (rrdetails), desiutil, desispec, redrock (all on PyPI).
  desitarget is optional (`[all]`): without it the page just lacks target-bit names.
- `MANIFEST.in` prunes `scratch/` (hundreds of MB of local test data) from sdists.

## Python API

**`prospect2.make_html(spectra, zcat=None, models=None, fits=None, ...)` is the entry point**
(`page.py`); the script only parses arguments, reads files and calls it. Goal (user, 2026-11-14):
users can read several coadds, filter and combine them, and make one page without temporary files.
User decisions:
- *Row-matched inputs:* `zcat`, `models` (`{band: [nspec, nwave_band]}`) and `fits` must match the
  spectra row for row (desispec convention, nothing silently reordered); TARGETIDs and shapes are
  checked, with a ValueError.
- *Spectra from `read_spectra(..., return_models=True, return_redshifts=True)` work as is:*
  `zcat`/`models` default to `spectra.redshifts`/`spectra.model`. Explicit arguments override them
  (logged).
- *A list of Spectra is accepted* (and lists of `zcat`/`models`/`fits` to match). `desispec.spectra.stack`
  keeps only the first input's `redshifts.meta`, so the TEMNAMnn/TEMVERnn template versions are
  compared *before* stacking, and differences raise (supporting mixed versions would need
  per-spectrum template bookkeeping; maybe later). A single pre-stacked Spectra is also accepted;
  then template consistency is the caller's responsibility.
- *`read_rrdetails` returns raw records*, one entry per zcat row (None where missing or not
  matching zcat), so lists from several files concatenate and filter like the zcats; the display
  lists, ≈-duplicate marks and sky exclusion are computed in `make_payload`.
- `make_payload` (dict) and `render_html` (str) are separate, so the payload can be tested without
  parsing HTML.
- *Only what's needed to plot is required* (user, 2026-11-15): wave/flux/ivar and a fibermap
  with TARGETID. TARGETID stays required because it row-matches the zcat and keys the VI, the
  filter and deep links; made-up IDs would let VI from different pages collide. Every other
  fibermap column is optional (`util.fibermap_value`): missing values become None/'' and show as
  "—"; `COADD_EXPTIME`/`COADD_FIBERSTATUS` fall back to the per-exposure `EXPTIME`/`FIBERSTATUS`;
  missing `COADD_NUMEXP` means 1 exposure (user); missing `OBJTYPE` means "target" if TARGETID ≥ 0;
  no `FIBER` shows the row as "Target <TARGETID>" / "#i", never as a fiber number; no coordinates
  hide the cutout and Legacy Viewer link.
- *The grouping decides the unique key* (`grouping.py`; user decisions 2026-11-15). DESI groups
  spectra by tile (SPGRP cumulative/pernight/perexp: a TARGETID can be on several tiles) or across
  tiles (healpix/uniqpix: one coadd per TARGETID per SURVEY/PROGRAM, never across them), and never
  mixes the two, so mixing raises. Before 2026-11-15 VI was keyed by TARGETID alone, which merged the
  VI of one target on two tiles.
  - Tile-based key: TILEID, FIBER, NIGHT (NIGHT = LASTNIGHT for cumulative; kept in the key even when
    unneeded, for bookkeeping, user choice), + EXPID when the fibermap has it (uncoadded spectra).
    Cross-tile key: SURVEY, PROGRAM, TARGETID. No SPGRP: tile-based if TILEID and FIBER are known,
    else TARGETID. The key is checked; if not unique, VI is keyed by row (`ROW`) with a warning.
  - Values come from fibermap columns, else headers. Headers are read *per input before stacking*:
    `stack` keeps at most the first input's headers. A single pre-stacked Spectra therefore loses
    the other inputs' NIGHT/SURVEY/PROGRAM (documented; pass a list).
  - PROGRAM: the primary header's (healpix grouping, e.g. `dark`), else the fibermap header's
    FAPRGRM (the tile's, e.g. `dark1b`). These differ: `main-dark-9901` (uniqpix, `dark`) includes
    tile 100614 (`DARK1B`). Harmless, since groupings are never mixed.
  - The VI export always has TARGETID, SURVEY, PROGRAM, then the remaining key columns (user: keep
    TARGETID on tile pages even though it isn't part of the key). It doesn't follow prospect1's VI
    format (user: not needed).
  - The JS builds `t.key` from `D.keycols` (nothing extra per target in the payload). Other spectra
    of the same TARGETID on the page are linked under the title ("also tile 2000 night …"); this is
    a Map built at load, no payload cost. Multi-tile pages label list rows `tile:fiber`.
- File names for the page come from `provenance=dict(source, redrock, model_source, rrdetails)`,
  all optional. They are shown without using screen space (user, 2026-11-15): as a tooltip on
  the page name in the header, and under "Made from" in the `?` help overlay (the header name is
  hidden on narrow screens).

## Tests

Four tiers in `py/prospect2/test/` (user-approved design, 2026-11-14):
1. *Unit tests* (no data): fiber lists, selection, output naming, q2 round trips, camera-model
   merging.
2. *Synthetic Spectra* (`desispec.test.util.get_blank_spectra` + the fibermap columns the page
   reads): `make_payload`/`make_html` behavior and errors, `parse()`. Models are passed in and
   alternate-fit models are off, so no redrock templates are needed.
3. *Real data* for specific fibers of petal 4, tile 1496 (masked cameras, NO DATA, alternate-fit
   rules, API page == script page), found via `$PROSPECT2_TEST_DATA`, `scratch/data/`, then the
   loa/nevis tile directory under `$DESI_ROOT`; skipped if absent. An explicitly set but wrong
   `$PROSPECT2_TEST_DATA` is an error, not a silent skip.
4. *Browser* (opt-in, `-m browser`): headless Chrome loads a synthetic page; no JS errors, and
   per-spectrum sums of the JS-decoded flux/model match the Python decoder.

Decisions:
- *No byte-identical golden page in pytest:* every intentional change would break it. `cmp`
  against the previous build stays the manual check for refactors (CLAUDE.md).
- *No real-data fixture in git:* 2–4 MB of binaries forever; tier 3 runs wherever the data are.
- `decode_q2` in `test/util.py` duplicates the JS `decodeQ2`; tier 4 guards against drift.
- Model round trip: error ≤ half a step. (Until 2026-11-14 the model was quantized with the exact
  scale but decoded with the 6-digit `mscale`, adding up to 5e-6·|model|; now `encode_model`
  quantizes with the rounded scale.)
- Chrome must run without `--user-data-dir`: with it, headless Chrome on macOS hangs.

## Data handling (Python)

**Cameras are merged with `desispec.coaddition.coadd_cameras`**, onto a single linear wavelength
grid. The JS assumes a linear grid (`wave0`, `dwave`), which is asserted in Python.

**The model is per camera, then merged with the flux's ivar weights** (`camera_models` /
`read_rrmodel` → `merge_camera_models`).
- prospect's `create_model` hard-cuts at mid-overlap, and the coadd_cameras model is an unweighted
  mean. Both give model≈0 (or half) where one camera is masked but the other has good data, e.g.
  fiber 2208 at 7537.6 Å.
- Where every camera is masked, the plain mean is used.

**rrmodel files are preferred when present.** `rrmodel-X.fits` next to `redrock-X.fits` is read,
with rows matched by TARGETID and the wavelength grids checked. A mismatch is an error (user
decision 2026-11-14; it used to fall back silently to computed models). Without an rrmodel file, or
with `--no-rrmodel`, models are computed from the templates + COEFF.
- *Known issue:* rrmodel models sag toward 0 within ~5 px of masked pixels. For example, fiber
  2183 z pix 809 is 22.1 in the rrmodel vs 29.6 computed, with data 29.4. This looks like masked
  pixels being zeroed before the resolution convolution.
- Elsewhere the two agree to float32 precision.

**Alternate redrock fits come from `rrdetails-X.h5`** (`read_rrdetails`, `--nfits N`, default 5,
user choice 2026-11-13). Each target's `zfit/<TARGETID>/zfit` table has 9 fits ranked by χ²; row 0
must match `redrock-X.fits` (z and SPECTYPE) or that target gets no alternates. The `zscan`
χ²(z) curves are not used yet.
- Stored compactly as `[z, spectype, subtype, χ² − best χ², zwarn, dup(, model row)]` arrays (~170 bytes per
  spectrum; JSON objects cost twice that and pushed the example page just under 10×).
- Shown Δχ² is relative to the *best fit*. Redrock's per-row `deltachi2` compares each fit with
  the next fit more than 1000 km/s away, so it is confusing when stepping through fits.
- Near-duplicates are kept (user choice), just faded and marked `≈#k`: `dup` = the first
  better-ranked fit within 1000 km/s (redrock's own criterion). E.g. a QSO fit 24 km/s from the best
  GALAXY fit has identical line markers but changes the key-line panels.
- **Alternate fits show their own redrock model when one is stored** (user choice 2026-11-14,
  replacing the earlier "hide the model for alternates"). Shifting the best-fit model to another
  fit's z would be misleading (e.g. a galaxy model slid to a QSO redshift), so the generator
  computes real models (`alt_fit_models`: rrdetails `coeff` + templates, each camera's resolution,
  merged like the best fit). Selection rule, measured to keep the 500-spectrum page at +3.5%:
  - only fits 2..N+1, `--alt-models N` (default 2);
  - only fits with χ² − best χ² < `--alt-max-dchi2` (default 100). Uncapped, 153 of the 500
    spectra had all 8 alternates under 100, almost all sky or ZWARN≠0 (noise fits);
  - never for OBJTYPE ≠ TGT (sky, unassigned/broken);
  - encoded like the best model but at `--alt-model-precision` (default σ/30, max error 0.017σ;
    staircases may show when zoomed + smoothed, extrapolated ~8–20 device px worst case vs 2.5–6
    at σ/100).
  Stored as payload `amodel`/`amscale` rows; a fit's array gets a 7th element = its row. Other
  alternates hide the model with the reason under the plot. The masked-model rules (MHIDE fill,
  no model in a fully masked camera) apply to these rows too (`AHIDE`).
  - Lower *wavelength* resolution was measured and rejected: the resolution-convolved model is
    barely sampled, so keeping every 2nd pixel gave >1σ errors at emission lines in 35% of
    spectra (to ~37σ) while saving only ~25%. Second-difference coding, coding as a residual from
    the best model, and shipping templates + coefficients to the browser were all worse.
- **Sky fibers get no alternate fits at all** (hard-coded, user choice 2026-11-14): redrock fits to
  sky residuals are noise.
- Selecting a fit sets the working z with `S.fit[i] = k`. `[ ] { }` nudges keep the fit; typing a
  z, `p` or ↺ return to the best fit. Key-line panels use the selected fit's spectype.
- ZWARN of each fit is only in the row tooltip: alternates mostly carry SMALL_DELTA_CHI2, so a
  warning icon on every row was clutter.

**The trial-z model is the pipeline model shifted in log λ**, not a refit. This ignores the change
in resolution, which is fine for display, and avoids shipping templates to the browser.

**Fiber lists use `10,20,30-40,50:60`.** Here `a-b` is inclusive and `a:b` is python-style
(excludes b). The syntax is the same in `-f` (`parse_fiberlist`) and in the page's filter box
(`parseIdList`).
- `parse_fiberlist` replaces desispec's `parse_fibers`, where `a-b` excludes b.
- `-f` and `-t` select the union.

**Output name defaults to `prospect-X.html` in the current directory** (from `coadd-X.fits` or
`spectra-X.fits`, `.gz` optional; `default_outfile`, `resolve_outfile`). User decisions,
2026-11-14:
- *Current directory, not next to the input:* production inputs are read-only.
- *Only for 1:1 pages:* the `prospect-X.html` name promises every spectrum of `coadd-X.fits`.
  With `-f`/`-t` (decided on the flags, not the result) or an input with another name, `-o FILE`
  is required, and the error says why. No auto-names like `prospect-X-f2000-2019.html`.
- *`-o` is cp-style:* an existing directory or a trailing `/` means a directory, which gets the
  default name; anything else is a file. An extension-based rule (no `.html` → directory) was
  rejected, because with auto-created directories a typo like `page.htlm` would make a directory.
  Instead, a file name not ending in `.html`/`.htm` gets a warning. `.html.gz` is not special:
  browsers don't decompress `file://` gz pages.
- *Missing directories are created* (DESI convention, cf. `desispec.io.util.makepath`, though most
  Unix tools would fail), with an `INFO: Created directory` line so typos are visible.
- *Existing files are overwritten* (regenerating is the common case), logged as `Overwrote
  existing …` instead of `Wrote …`.
- Naming errors are raised before `read_spectra`, so they are immediate.

## Payload contents

**The payload is internal** (user decision, 2026-11-15): one JSON object, built by
`make_payload()` (`page.py`), written by `render_html()` in place of `"__DATA__"`, and read by the
viewer as `const D`. It always ships in the same file as the JS that reads it, so it has no
version contract and can change freely; keep this section current instead.

Top level (`D.x`):

| key | meaning |
|---|---|
| `title` | page title: `title` argument, else `provenance['source']`, else `'prospect2'` |
| `source` | input file name, else `title`; never null. Keys the VI in localStorage and names the VI export |
| `redrock`, `model_source`, `rrdetails` | provenance file names (or `'input models'`, `'redrock templates + COEFF'`), or null; shown in the header tooltip and help overlay |
| `wave0`, `dwave`, `nwave` | the merged linear wavelength grid (`coadd_cameras`); JS `WAVE`, `NW` |
| `nspec`, `flux_unit` | number of spectra; flux unit label |
| `spgrp`, `keycols` | the inputs' SPGRP (or null) and the columns of the unique key that keys the VI (`grouping.key_columns`) |
| `targets` | one dict per spectrum, row order = data rows (below) |
| `lines` | emission/absorption line lists (`io.read_lines`) |
| `overlaps` | camera overlap regions `[lo, hi, label]` |
| `enc`, `noise`, `flux`, `model`, `mscale` | the q2 arrays (next section); `model`/`mscale` null without a zcat |
| `amodel`, `amscale`, `alt_models` | only when alternate-fit models exist: their rows, and the `--alt-models`/`--alt-max-dchi2` settings |

`targets[i]` (JS `T[i]`), all plain JSON:
- Always: `i`, `targetid` (**string**, int64 > 2^53), `fiber`, `tileid`, `night`, `expid`, `survey`,
  `program` (grouping values, `grouping.grouping_columns`), `ra`, `dec`, `objtype`,
  `morph`, `classes`/`primary` (target-bit names), `mag` (`{g, r, z, w1, w2}`), `fibermag_r`,
  `ebv`, `exptime`, `numexp`, `fiberstatus`, `fiberstatus_names`, `snr` (`{band: median S/N}`),
  `masked_cams` (`[lo, hi, band]` for cameras with no unmasked pixels). Values from optional
  fibermap columns can be null (strings: `''`), so the JS must handle that; see "Python API".
- With a zcat: `z`, `zerr`, `zwarn`, `zwarn_names`, `spectype`, `subtype`, `deltachi2`, `chi2`,
  `npixels`. Absent, not null, without one. Redrock's no-fit placeholders are kept raw; see
  "Degenerate spectra".
- With rrdetails, non-sky only: `fits`, compact arrays `[z, spectype, subtype, Δχ², zwarn, dup(,
  model row)]` (`_fit_table`). The JS expands them to objects `{z, st, sub, dchi2, zw, dup, m}`
  at load, and adds `t.ngood` (unmasked pixel count) after decoding.

Decoded arrays: `FLUX`, `NOISE`, `MODEL` are `Float32Array(nspec·NW)`, row *i* = `T[i]`
(`row(arr, i)`), masked = NaN. `AMODEL` rows are indexed by a fit's `m`. `MHIDE`/`AHIDE` flag the
model pixels to hide (see "Masked model rule").

**Recipes:**
- *New per-target value:* add it to the `t = dict(...)` in `make_payload`, reading fibermap
  columns with `col(...)` (`fibermap_value`) so a missing column gives null, and read `T[i].key` in
  the JS; no template placeholder is needed. Show it in `renderInfo()`; to filter on it, add it to
  `FIELDS` (numeric comparisons) or to the label string in `compileTerm` (substring match), and
  update the filter help in the template and README. Test it in `test_page.py`. Nothing checks
  that a payload key is used by the JS, or that a key the JS reads exists (it is just
  `undefined`).
- *New per-pixel array:* encode it in `encode.py`, decode it in `decodeQ2()` **and** in
  `test/util.py` `decode_q2`, then run `pytest -m browser`.

## Payload encoding ("q2")

The payload is gzip + base64, decoded in the browser with `DecompressionStream`. The example page
is 485 kB, 10.2× smaller than the 4.95 MB Bokeh page for the same 20 spectra.

- **Packing (`_esc8`):** every integer array is stored as an int8 stream plus an int32 escape
  stream (-128 = take the next escape value; for flux, -127 = masked → NaN). gzip compresses
  bytes much better than int16/int32: this was ~20% smaller than the earlier int16/int32 "q1"
  packing of the *same* quantized values (2026-11-13; decoded floats verified bit-identical), and
  decoding is slightly faster (77 vs 85 ms for 500 spectra).
- **Noise** is stored as `round(100·ln σ)`, delta-coded, i.e. 1% steps.
- **Flux** is quantized at 0.1σ against the *decoded* noise, so the error is ≤ 0.05σ everywhere.
- **Model** is delta-coded, in steps of (median *data* σ)/N, with `--model-precision N`
  defaulting to 100. If a spectrum has no unmasked data, the step is span(model)/(10N).
  - The model is noiseless, so coarse steps show as staircases. σ/10 was visibly stepped.
  - At N=100, steps are invisible in the full view but can reach 3–6 device px when zoomed and
    heavily smoothed.
  - N=1000 keeps steps < 0.4 px, at +1.2 MB per 500 spectra. The user chose 100 as the default.
  - A model-only rule (span/N) would need span/5e4 for the same worst case, because span/σ varies
    from 0.7 to 118.
- *Further size options measured (2026-11-13), not done:* gzip-packing the app JS + uPlot inside
  the page (saves ~60–70 kB fixed per page, +1–3 ms load; only matters for small pages; needs
  global-scope `<script>` injection, `//# sourceURL`, and a `--no-pack` flag to keep testing and
  debugging workable); noise in 2% steps (−1.5 kB/spectrum, lossy but invisible); flux at 0.2σ
  (−1.2 kB/spectrum); model precision 30 (−0.8 kB/spectrum, but reintroduces visible staircases,
  rejected); base85 (~−6%); brotli/xz (~−17%, needs an inlined JS decoder).
- **TARGETIDs are strings.** int64 values above 2^53 lose precision as JSON numbers (…322 became
  …320). Negative TARGETIDs are sky or unassigned fibers.

## Plotting (JS)

**Y-range is robust**: the 1–99th percentile of the flux in the visible window, padded, then
widened to include the model's full min/max (not percentiles) and 0. A "Full" toggle uses the flux's full range instead of the percentiles.
- Masked data never enters the y-range.
- uPlot gotcha: use `y: {auto: true, range: u => calcY(u)}`. `calcY` reads
  `u.series[0].min/max`, because the range function runs before `scales.x` updates. For refits
  without an x change, set `series[0].min/max` and call `setScale("y", …)` (`fitY`).

**Smoothing**: flux smoothing is ivar-weighted Gaussian. Model smoothing is plain Gaussian, by user
choice.

**Masked model rule** (user-specified):
- Only where the data are masked (ivar=0) *and* the pipeline model == 0, the model pixel is
  flagged.
- It is filled by linear interpolation (`fillMaskedModel`) so that shifting and smoothing don't
  see fake zeros, then hidden after shift and smooth.
- Masked pixels with model ≠ 0 are neither interpolated nor hidden.

**Camera overlaps are shaded** (computed from the per-camera grids, key `c`), instead of sky or
telluric bands, per user preference.

**Line markers** are SDSS-style short ticks (`v` cycles Full | Data | Model | Edge; default Data).
Emission ticks sit above the spectrum and absorption ticks below.
- *Tick placement:* the tick is anchored outside the extreme of *everything plotted* near the line.
  The noise trace counts for emission only, because it sits near 0 and would drag absorption
  labels to the bottom.
- *Label placement:* the label is pushed clear of the data under its full width, otherwise
  unsmoothed noise spikes sit on labels. Labels are also clamped inside the plot horizontally.
- *Merging:* same-name lines within 60 px share a bracketed label (doublets). Any lines within
  10 px merge into one label ("K/H").
- *Spacing:* there is an 8 px gap from the data (`&lg=`), and ticks carry a surface-coloured halo.

**Marker colours** default to ink: emission ink, absorption dark gray. Other palettes are plum,
crimson and violet (`V`, `&lc=`).
- The user found the original violet too close to the blue flux.
- Measured ΔE (OKLab×100) under colour-vision-deficiency simulation: violet vs the blue flux is ≈5
  (deutan), and any green vs the orange model is ≈1 (protan). Neutral ink is ≥25 from both in
  every simulation.
- Data colours are flux `#2a78d6` and model `#eb6834`.

**Text goes in uPlot's `draw` ("over") layer.** `drawClear` ("under") runs before the grid, so
gridlines strike through text drawn there.

## Degenerate spectra

**Fully masked spectra** (every pixel ivar=0, COADD_NUMEXP=0; 6 in petal 4, e.g. fiber 2036):
- They get a NO DATA badge and a hatched plot with a reason box (NUMEXP, FIBERSTATUS). There are
  no key-line panels and the grid thumbnail is hatched.
- An empty plot otherwise looks like a rendering bug.
- Redrock's placeholders (NPIXELS=0, CHI2=9e99, DELTACHI2≈1.9e84) are shown as "no fit" / "—" by
  `noFit()`. They are kept raw in the Redrock table and never sorted or filtered on as real values.
- No model is drawn for these spectra.

**OBJTYPE '' + negative TARGETID → UNASSIGNED**, not SKY. OBJTYPE SKY can also have a negative
TARGETID (stuck sky fibers).

**An entirely masked camera** (e.g. BADAMPZ, fibers 2026–2028 z, 2277 b+z) gets special treatment:
- The coadd stores flux=ivar=0 there, so there are no masked measurements to show.
- The span is hatched and labelled "z camera masked (BADAMPZ)", with the reason taken from
  FIBERSTATUS.
- Line markers are suppressed inside the span, and label boxes are reserved so that line labels
  avoid it.
- **The model is never drawn there.** The user tried a light "ghost" model and rejected it: a PCA
  model extrapolated over a whole camera is distractingly wrong.
- Small masked regions remain plain gaps.

## UI

- There are inspect and grid views. Key-line zoom panels are chosen by SPECTYPE; clicking a panel
  jumps there.
- The right sidebar has no separate "Redrock" key/value section (dropped 2026-11-13, user choice, to
  save vertical space): it duplicated the header facts. What was unique moved to the header: ZWARN
  bit names are spelled out in muted text below the Quality badge (absolutely positioned in a
  space the facts always reserve, so the layout doesn't jump when ZWARN changes; user request), and Χ² / Npix is the ΔΧ² fact's tooltip. If the
  bit names ever crowd the header, the badge tooltip alone is acceptable to the user. The "Redrock
  fits" table (rrdetails) takes the section's place.
- The legend is replaced by a hover readout under the plot. Both sidebars are collapsible (`b`/`i`)
  and the choice is remembered in localStorage.
- Filter syntax: `z>1.5 0.79<z<0.81 dchi2<25 ELG !QSO nodata unassigned 2030-2040`.
  Whitespace-separated terms must all match (AND); `!` negates a term. By user choice (2026-11-13),
  terms never contain spaces: `z > 1.5` is not supported and fiber lists must be written without
  spaces (`10,20`, not `10, 20`), so tokenizing is a plain whitespace split. Terms that contain
  `<`, `>` or `=` but don't parse as a comparison are listed in red under the box and ignored
  (not treated as "match nothing"). Commas only join fiber/TARGETID lists, not terms.
- Deep links: `#i=7&s=3&view=grid&theme=dark&q=<filter>&ls=&lc=&lg=`. The hash is kept current
  on select and filter.
- VI (visual inspection, quality 0–4) is stored in localStorage per source file, each entry keyed
  by the spectrum's unique key (see "Python API"), and exported as CSV. There is no server.
- Sort by fiber sorts by tile, then fiber. The filter has `tile` and `night` fields, and survey and
  program are text-matched.
- Light and dark themes are supported. The canvas reads CSS custom properties, so re-render on
  theme change (`rethemed`). The ◐ button flips to the opposite of the theme currently shown (pages start
  by following the OS setting). It used to cycle auto → light → dark, which needed two clicks to
  change anything visible when the OS was in light mode (user report, Safari).

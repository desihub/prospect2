# prospect2 implementation notes

Key design decisions and the reasons behind them. `README.md` covers usage; this file covers *why*.
Keep it current when a decision changes, and record the reason (often a user preference or a
specific fiber that exposed a problem).

## Architecture

**Python builds a data payload; a static hand-written JS app renders it.** `make_prospect2.py`
reads the coadd/redrock files and writes one self-contained HTML file: `prospect2_template.html`
with uPlot 1.6.31 (`vendor/`, ~50 KB) and the payload inlined.

- *Why not patch Bokeh (the original prospect)?* The Bokeh page was 4.7 MB for 20 spectra plus
  CDN-loaded BokehJS, and its interaction model is hard to customize.
- *Why not Plotly?* The library alone is ~3.5 MB.
- *Why not Panel/Dash?* They need a running server; pages must work offline as a file, by email,
  or on a static web server.
- *Why uPlot?* It is tiny and fast for 500 × 7781-pixel spectra, and exposes canvas draw hooks for
  custom overlays (line markers, hatching).

**No dependency on prospect.** The line lists (`data/{emission,absorption}_lines.txt`) and
`load_redrock_templates` were copied from prospect, so the prototype stands alone. The
dependencies are numpy, astropy, desispec, redrock (+ templates via `$RR_TEMPLATE_DIR`), and
optionally desitarget. License: BSD 3-clause, same authors as prospect.

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
with rows matched by TARGETID and the wavelength grids checked. On any mismatch, or with
`--no-rrmodel`, models are computed from the templates + COEFF.
- *Known issue:* rrmodel models sag toward 0 within ~5 px of masked pixels. For example, fiber
  2183 z pix 809 is 22.1 in the rrmodel vs 29.6 computed, with data 29.4. This looks like masked
  pixels being zeroed before the resolution convolution.
- Elsewhere the two agree to float32 precision.

**Alternate redrock fits come from `rrdetails-X.h5`** (`read_rrdetails`, `--nfits N`, default 5,
user choice 2026-11-13). Each target's `zfit/<TARGETID>/zfit` table has 9 fits ranked by χ²; row 0
must match `redrock-X.fits` (z and SPECTYPE) or that target gets no alternates. The `zscan`
χ²(z) curves are not used yet.
- Stored compactly as `[z, spectype, subtype, χ² − best χ², zwarn, dup]` arrays (~170 bytes per
  spectrum; JSON objects cost twice that and pushed the example page just under 10×).
- Shown Δχ² is relative to the *best fit*. Redrock's per-row `deltachi2` compares each fit with
  the next fit more than 1000 km/s away, so it is confusing when stepping through fits.
- Near-duplicates are kept (user choice), just faded and marked `≈#k`: `dup` = the first
  better-ranked fit within 1000 km/s (redrock's own criterion). E.g. a QSO fit 24 km/s from the best
  GALAXY fit has identical line markers but changes the key-line panels.
- **The model is hidden for alternate fits** (user choice "a"), with a note under the plot. Only
  the best-fit model is stored; shifting it to another fit's z would be misleading, e.g. a galaxy
  model slid to a QSO redshift. Storing real models for alternates is possible (~+2 kB/spectrum
  each, the generator can compute them from templates + COEFF) but was not wanted.
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

**Y-range is robust**: the 1–99th percentile of flux and model in the visible window, plus
padding, including 0. A "Full" toggle switches to the full range.
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
- VI (visual inspection, quality 0–4) is stored in localStorage keyed by source file and exported
  as CSV. There is no server.
- Light and dark themes are supported. The canvas reads CSS custom properties, so re-render on
  theme change (`rethemed`). The ◐ button flips to the opposite of the theme currently shown (pages start
  by following the OS setting). It used to cycle auto → light → dark, which needed two clicks to
  change anything visible when the OS was in light mode (user report, Safari).

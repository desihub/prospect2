# prospect2 prototype

```
python make_prospect2.py -i coadd-4-1496-thru20240214.fits -o prospect-4-1496-thru20240214.html
```

Matching `redrock-*.fits`, `rrmodel-*.fits` and `rrdetails-*.h5` files are found next to the coadd.
Use `-f`/`-t` to select a subset (see below).

Page sizes for petal 4, tile 1496:

| spectra | prospect (Bokeh)                 | prospect2           |
|---------|----------------------------------|---------------------|
| 20      | 4.95 MB (+ Bokeh JS from CDN)    | 0.51 MB, offline    |
| 500     | 107 MB (+ Bokeh JS from CDN)     | 8.87 MB, offline    |

Deep links: `#i=7&s=3` (target index, smoothing), `&view=grid`, `&theme=dark`, `&q=<filter>`.

Fiber/TARGETID lists, same syntax in `-f` and in the page's filter box: `10,20,30-40,50:60` =
fibers 10, 20, 30–40 (inclusive), 50–59 (python-style). In the filter box, integers ≥ 5000 or
negative are TARGETIDs (exact match, or partial match for a single number). `-f` and `-t` together
select the union.

Filter box: whitespace-separated terms, all of which must match (`2030-2040 ELG 0.79<z<0.81 !QSO`).
A term is a comparison on `z`, `dchi2`, `fiber`, `snr` or `zwarn` (`z>1.2`, `1.2<=z`, `0.79<z<0.81`),
a fiber/TARGETID list (no spaces inside it), or a case-insensitive text match on the target's labels
(`ELG`, `GALAXY`, `nodata`, `unassigned`). Prefix `!` negates a term. Terms must not contain spaces
(`z > 1.2` doesn't work). Terms that can't be parsed are shown in red under the box and ignored.

Dependencies: numpy, astropy, desispec, redrock (+ templates via `$RR_TEMPLATE_DIR`), and
optionally desitarget for target-bit names. No dependency on prospect: the line lists
(`data/{emission,absorption}_lines.txt`) and `load_redrock_templates` were copied from it.
License: BSD 3-clause (`LICENSE.rst`), same as prospect.

## Architecture

* Python (`make_prospect2.py`) builds a JSON payload; the HTML template (`prospect2_template.html`)
  is a small hand-written JS app plus inlined uPlot (50 KB, `vendor/`). No Bokeh, no server.
* Cameras merged with `desispec.coaddition.coadd_cameras`. Per-camera redrock models come from the
  matching `rrmodel*.fits` (`redrock-X.fits` → `rrmodel-X.fits`, rows matched by TARGETID, wavelength
  grids checked) if present, else are computed from templates + COEFF (`camera_models`; forced with
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
* VI stored in localStorage, exported as CSV.
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

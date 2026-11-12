# prospect2 prototype

```
python alt/make_prospect2.py -i data/coadd-4-1496-thru20240214.fits -o alt/prospect2_example.html -f 2000:2020
```

| file | spectra | size |
|---|---|---|
| `../example_prospect.html` (current Bokeh prospect) | 20 | 4.7 MB (+ Bokeh JS from CDN) |
| `prospect2_example.html` | 20 | 0.55 MB, fully offline |
| `prospect2_petal4_all500.html` | 500 | 10.6 MB, fully offline (12.2 MB with `--model-precision 1000`) |

Deep links: `#i=7&s=3` (target index, smoothing), `&view=grid`, `&theme=dark`, `&q=<filter>`.

Fiber/TARGETID lists, same syntax in `-f` and in the page's filter box: `10,20,30-40,50:60` =
fibers 10, 20, 30–40 (inclusive), 50–59 (python-style). In the filter box, integers ≥ 5000 or
negative are TARGETIDs (exact match, or partial match for a single number); terms separated by
spaces must all match (`2030-2040 ELG z>1.2`). `-f` and `-t` together select the union.

Dependencies: numpy, astropy, desispec, redrock (+ templates via `$RR_TEMPLATE_DIR`), and
optionally desitarget for target-bit names. No dependency on prospect: the line lists
(`data/{emission,absorption}_lines.txt`) and `load_redrock_templates` were copied from it.
License: BSD 3-clause (`LICENSE.rst`), same as prospect.

## Architecture

* Python (`make_prospect2.py`) builds a JSON payload; the HTML template (`prospect2_template.html`)
  is a small hand-written JS app plus inlined uPlot (50 KB, `vendor/`). No Bokeh, no server.
* Cameras merged with `desispec.coaddition.coadd_cameras`. The redrock model is computed per camera
  (`camera_models`) and merged with the same ivar weights as the flux (`merge_camera_models`).
  prospect's `create_model` hard-cuts at mid-overlap, giving model=0 where that camera is masked
  but the other camera has good data (e.g. fiber 2208 at 7537.6 Å).
* **Data encoding ("q1")**: noise as 1%-step log, delta-coded int16; flux quantized at 0.1σ
  (int16, int32 only if S/N > 3276/pix); model delta-coded int32 in steps of (median data σ)/N,
  `--model-precision N` (default 100; span(model)/(10N) for spectra with no unmasked data).
  The model is noiseless, so steps can show as staircases: N=100 is invisible in the full view but
  can reach ~3–6 device px in zoomed, heavily smoothed flat regions; N=1000 stays < 0.4 px
  (+1.6 MB per 500 spectra). Flux error ≤ 0.05σ, mask preserved.
  ~8.5x smaller than the Bokeh page for the same 20 spectra.
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

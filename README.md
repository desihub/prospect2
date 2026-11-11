# prospect2 prototype

```
python alt/make_prospect2.py -i data/coadd-4-1496-thru20240214.fits -o alt/prospect2_example.html -f 2000-2020
```

| file | spectra | size |
|---|---|---|
| `../example_prospect.html` (current Bokeh prospect) | 20 | 4.7 MB (+ Bokeh JS from CDN) |
| `prospect2_example.html` | 20 | 0.55 MB, fully offline |
| `prospect2_petal4_all500.html` | 500 | 10.6 MB, fully offline (12.2 MB with `--model-precision 1000`) |

Deep links: `#i=7&s=3` (target index, smoothing), `&view=grid`, `&theme=dark`.

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

# prospect2 (alt/)

Prototype replacement for the DESI prospect spectral viewer. A Python generator writes a single
offline HTML page from a hand-written JS app and inlined uPlot.

- **Read `docs/implementation.md` before changing behavior.** It records the design decisions and
  the reasons for them, many of which were explicit user choices. Don't silently reverse one. If a
  decision changes, update that file with the new reason.
- `README.md` covers usage, flags and the deep-link syntax. Keep it in sync with user-visible
  changes.

## Layout

- `make_prospect2.py`: the generator. It reads coadd + redrock (+ optional rrmodel), merges the
  cameras, encodes the payload ("q2"), and fills the template.
- `prospect2_template.html`: the whole viewer (CSS + JS, ~1200 lines). The placeholders
  `/*__UPLOT_CSS__*/`, `/*__UPLOT_JS__*/`, `"__DATA__"` and `__TITLE__` are filled by the
  generator.
- `vendor/`: uPlot 1.6.31 (MIT, see `uPlot.LICENSE`). Don't upgrade casually; the y-range hooks
  depend on its behavior.
- `data/`: emission/absorption line lists (copied from prospect).
- Example inputs live in `../data/` (petal 4, tile 1496: coadd, redrock, rrmodel). They are not in
  this repo.

## Running

Use the `desi` conda env's python directly. `conda run` drops heredoc stdin.

```
export RR_TEMPLATE_DIR=~/desi/git/redrock-templates
cd ~/desi/dev/prospect2
~/miniforge3/envs/desi/bin/python alt/make_prospect2.py -i data/coadd-4-1496-thru20240214.fits \
    -o alt/prospect2_example.html -f 2000:2020
~/miniforge3/envs/desi/bin/python alt/make_prospect2.py -i data/coadd-4-1496-thru20240214.fits \
    -o alt/prospect2_petal4_all500.html
```

Note that `-f 2000-2020` is *inclusive* (21 fibers); this differs from desispec. After changing
the generator or the template, regenerate both example pages.

## Testing

There is no test suite. Verify by rendering:

- Screenshot with headless Chrome:
  `"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --window-size=1600,1100 --virtual-time-budget=8000 --screenshot=OUT.png "file://…/prospect2_example.html#i=7"`
- Add `--enable-logging=stderr --v=0` and grep for `CONSOLE` to catch JS errors and `console.log`
  output.
- For interaction tests, inject a `<script>` into a *scratch copy* of the page that dispatches
  keydown events or sets state. Do the injection with Python string replacement, not sed: `#` and
  `/` in the JS break sed.
- Look at the screenshots: check label collisions, clipping and overlap with the data. Useful
  fibers in petal 4:

  | fiber(s) | what it exercises |
  |---|---|
  | 2036 | fully masked / NO DATA |
  | 2026–2028 | z camera masked |
  | 2277 | b and z cameras masked |
  | 2208 | camera masked in an overlap |
  | 2183 | rrmodel sag next to a masked gap |

- When the output should not change, compare the generated HTML byte-for-byte (`cmp`) against
  the previous build.

## Conventions and pitfalls

- TARGETIDs must go to JS as **strings** (int64 > 2^53).
- Never draw the model where the whole camera is masked, nor for spectra with no redrock fit
  (`noFit`).
- Masked pixels must never enter the y-range.
- Draw text in uPlot's "over" layer, not in `drawClear`.
- Colours: data are blue (flux) and orange (model). Check new marker or overlay colours for
  colour-vision-deficiency separation (ΔE) against both, not by eye. Support light and dark themes
  (CSS custom properties; canvas code reads them via `css()`).
- `[hidden]` needs `display:none !important`.
- Redrock SUBTYPE can be `--` (masked); strip it.

## Git

- Commit only when asked.
- Generated `prospect2_*.html` pages are git-ignored (except `prospect2_template.html`); don't add
  them.

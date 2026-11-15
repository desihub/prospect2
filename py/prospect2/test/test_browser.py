#
# See top-level LICENSE.rst file for Copyright information

"""
Tier 4 (opt-in: pytest -m browser): load a page in headless Chrome, check that
there are no JS errors, that the browser's decodeQ2() gives the same flux and
model as the Python decoder, and that the filter box's parseIdList() reads
fiber lists like util.parse_fiberlist (-f).  Chrome comes from $CHROME, the macOS app, or
google-chrome / chromium on $PATH; skipped if none is found.
"""

import os, re, json, shutil, subprocess
import numpy as np
import pytest

from prospect2.page import make_html
from prospect2.util import parse_fiberlist
from prospect2.test.util import synthetic_spectra, synthetic_zcat, synthetic_models, decode_q2, payload_from_html

pytestmark = pytest.mark.browser

#- fiber lists that the filter box and -f must read the same way
FIBERLISTS = ['10', '10,20', '30-40', '50:60', '10,20,30-40,50:60', '0', '4999', '7,', '1-1', '5:6', '3,3']
#- ... and single numbers that the filter box (only) reads as TARGETIDs: >= 5000 or < 0
TIDTOKENS = ['5000', '39628266889283065', '-1']

#- injected after the app: waits for decoding, then logs per-spectrum sums and parsed fiber lists
CHECK_JS = """<script>
const P2LISTS = __FIBERLISTS__, P2TIDS = __TIDTOKENS__;
setTimeout(() => {
  const sums = a => Array.from({length: D.nspec}, (_, s) => {
    let sm = 0, nnan = 0;
    for (let k = 0; k < NW; k++) { const v = a[s * NW + k]; if (Number.isFinite(v)) sm += v; else nnan++; }
    return [sm, nnan];
  });
  const ids = tok => { const r = parseIdList(tok); return r && { fibers: r.fibers, tids: [...r.tids] }; };
  console.log("P2CHECK " + JSON.stringify({flux: sums(FLUX), model: sums(MODEL),
                                           lists: P2LISTS.map(ids), tids: P2TIDS.map(ids)}));
}, 3000);
</script>
</body>"""


def _chrome():
    for c in (os.environ.get('CHROME'), '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
              shutil.which('google-chrome'), shutil.which('chromium'), shutil.which('chromium-browser')):
        if c and os.path.exists(c):
            return c
    return None


def test_page_loads_and_decodes(tmp_path):
    chrome = _chrome()
    if chrome is None:
        pytest.skip('no Chrome found (set $CHROME)')
    sp = synthetic_spectra(nspec=3)
    sp.ivar['z'][1] = 0                   #- a masked camera
    for b in sp.bands:                    #- and a spectrum with no data
        sp.ivar[b][2] = 0
    html = make_html(sp, zcat=synthetic_zcat(sp), models=synthetic_models(sp, 1.5), alt_models=0)
    page = tmp_path / 'page.html'
    check = CHECK_JS.replace('__FIBERLISTS__', json.dumps(FIBERLISTS)).replace('__TIDTOKENS__', json.dumps(TIDTOKENS))
    page.write_text(html.replace('</body>', check, 1), encoding='utf-8')

    #- no --user-data-dir: with it, headless Chrome on macOS hangs instead of exiting
    cmd = [chrome, '--headless=new', '--disable-gpu',
           '--window-size=1600,1100', '--virtual-time-budget=8000', '--enable-logging=stderr', '--v=0',
           f'--screenshot={tmp_path / "page.png"}', f'file://{page}#i=1']
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stderr
    console = [line for line in out.splitlines() if 'CONSOLE' in line]
    errors = [line for line in console if 'Uncaught' in line or 'Error' in line]
    assert not errors, '\n'.join(errors)
    m = re.search(r'"P2CHECK (.*)", source', out)
    assert m, 'page did not report decoded data:\n' + '\n'.join(console)
    js = json.loads(m.group(1))

    flux, _, model, _ = decode_q2(payload_from_html(html))
    for name, arr in (('flux', flux), ('model', model)):
        a = arr.astype(np.float32)            #- JS decodes into Float32Arrays
        expected = [[float(np.nansum(r.astype(np.float64))), int(np.isnan(r).sum())] for r in a]
        for (jsum, jnan), (psum, pnan) in zip(js[name], expected):
            assert jnan == pnan, f'{name}: masked pixel counts differ'
            assert np.isclose(jsum, psum, rtol=1e-6, atol=1e-6), f'{name}: sums differ'

    #- fiber lists: JS ranges expand to the same fibers as parse_fiberlist
    for spec, ids in zip(FIBERLISTS, js['lists']):
        assert ids is not None and not ids['tids'], spec
        jsfibers = sorted(f for a, b in ids['fibers'] for f in range(a, b + 1))
        assert jsfibers == sorted(parse_fiberlist(spec).tolist()), spec
    #- documented difference: in the filter box, numbers >= 5000 or < 0 are TARGETIDs
    for tok, ids in zip(TIDTOKENS, js['tids']):
        assert ids == {'fibers': [], 'tids': [tok]}, tok

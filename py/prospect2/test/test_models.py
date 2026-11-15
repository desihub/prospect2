#
# See top-level LICENSE.rst file for Copyright information

"""Tier 1: merging per-camera models with the flux's ivar weights"""

import numpy as np

from desispec.coaddition import coadd_cameras
from prospect2.models import merge_camera_models
from prospect2.test.util import synthetic_spectra


def test_merge_camera_models():
    sp = synthetic_spectra(nspec=2)
    wave = coadd_cameras(sp).wave['brz']
    models = {b: np.full(sp.flux[b].shape, v) for b, v in zip('brz', (1.0, 3.0, 5.0))}
    #- b/r overlap: spectrum 0 weights r 3x more; spectrum 1 has r masked at the overlap pixel
    lam = sp.wave['r'][0]
    jb = np.argmin(np.abs(sp.wave['b'] - lam))
    sp.ivar['r'][0, 0] = 12.0
    sp.ivar['r'][1, 0] = 0.0
    #- z camera fully masked for spectrum 1 beyond the r/z overlap: plain mean only where all masked
    sp.ivar['b'][1, jb] = 0.0
    m = merge_camera_models(sp, models, wave)
    j = np.argmin(np.abs(wave - lam))
    assert np.isclose(m[0, j], (4 * 1 + 12 * 3) / 16)       # ivar-weighted
    assert np.isclose(m[1, j], 2.0)                          # both masked -> plain mean of 1 and 3
    assert np.isclose(m[0, 0], 1.0) and np.isclose(m[0, -1], 5.0)

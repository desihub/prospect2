#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.models
================

Redrock models: computed from templates + coefficients, convolved with each
camera's resolution, and merged across cameras with the flux's ivar weights.
"""

import os
import numpy as np


def load_redrock_templates(template_dir=None, zcat_header=None):
    """
    Load redrock templates as a dict keyed by (SPECTYPE, SUBTYPE), silencing
    redrock's stdout.  If zcat_header has TEMNAMnn/TEMVERnn keywords, load the
    template versions used for the fit.  (From prospect.utilities.)
    """
    import contextlib
    import redrock.templates
    with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
        if zcat_header is not None:
            return redrock.templates.load_templates_from_header(
                zcat_header, template_dir=template_dir, asdict=True)
        return redrock.templates.load_templates(template_dir, asdict=True)


def camera_models(sp, zcat):
    """
    Redrock best-fit model per camera, convolved with each camera's resolution
    (same recipe as prospect.viewer.create_model, but without its hard cut at
    the middle of each camera overlap).  Returns {band: [nspec, nwave_band]}.
    """
    from numpy.ma.core import MaskedConstant
    from desispec.interpolation import resample_flux
    hdr = zcat.meta if 'TEMNAM00' in zcat.meta else None
    templates = load_redrock_templates(zcat_header=hdr)
    out = {b: np.zeros(sp.flux[b].shape) for b in sp.bands}
    for i, zb in enumerate(zcat):
        subtype = zb['SUBTYPE']
        if isinstance(subtype, MaskedConstant):
            subtype = zcat['SUBTYPE'].data.data[i]
        if isinstance(subtype, bytes):
            subtype = subtype.decode('utf-8')
        tx = templates[(zb['SPECTYPE'], subtype)]
        restmodel = tx.flux.T.dot(zb['COEFF'][:tx.nbasis]).T
        for b in sp.bands:
            mx = resample_flux(sp.wave[b], tx.wave * (1 + zb['Z']), restmodel)
            out[b][i] = sp.R[b][i].dot(mx)
    return out


def merge_camera_models(sp, models, wave):
    """
    Merge per-camera models onto the coadd_cameras grid with the same ivar weights
    used for the flux, so a camera that is masked at a pixel (model==0 there) does
    not pull the model down.  Where every camera is masked, use the plain mean.
    """
    dw = wave[1] - wave[0]
    num = np.zeros((sp.num_spectra(), len(wave)))
    den = np.zeros_like(num); msum = np.zeros_like(num); cnt = np.zeros_like(num)
    for b in sp.bands:
        j0 = int(round((sp.wave[b][0] - wave[0]) / dw))
        sl = slice(j0, j0 + len(sp.wave[b]))
        assert np.allclose(wave[sl], sp.wave[b]), f'camera {b} not aligned with merged grid'
        iv = sp.ivar[b] * ((sp.mask[b] == 0) if sp.mask is not None else 1)
        num[:, sl] += iv * models[b]; den[:, sl] += iv
        msum[:, sl] += models[b]; cnt[:, sl] += 1
    return np.where(den > 0, num / np.where(den > 0, den, 1), msum / np.maximum(cnt, 1))


def alt_fit_models(sp, wave, zcat, fits, fitlists, nalt=2, max_dchi2=100):
    """
    Redrock models for alternate fits, computed from the rrdetails templates + coeff, convolved with
    each camera's resolution and merged like the best-fit model.  Only fits ranked 2..nalt+1 with
    chi2 - best chi2 < max_dchi2 get a model; sky and unassigned fibers (OBJTYPE != TGT) never do.
    `fits` is the per-spectrum list of raw rrdetails records (see prospect2.io.read_rrdetails) and
    `fitlists` the matching display lists; appends the model's row index to each such fit's entry
    in `fitlists` (in place).
    Returns (models [nrow, nwave], spectrum index of each row), or (None, None) if there are none.
    """
    from desispec.interpolation import resample_flux
    rows = []
    for i, d in enumerate(fits):
        if d is None or fitlists[i] is None:
            continue
        if str(sp.fibermap['OBJTYPE'][i]).strip() != 'TGT':
            continue
        for k in range(1, min(nalt + 1, len(d))):
            if d['chi2'][k] - d['chi2'][0] < max_dchi2:
                fitlists[i][k].append(len(rows))
                rows.append((i, d[k]))
    if not rows:
        return None, None
    hdr = zcat.meta if 'TEMNAM00' in zcat.meta else None
    templates = load_redrock_templates(zcat_header=hdr)
    idx = np.array([i for i, _ in rows])
    sub = sp[idx]
    cam = {b: np.zeros(sub.flux[b].shape) for b in sp.bands}
    s = lambda x: (x.decode() if isinstance(x, bytes) else str(x)).strip()
    for j, (i, d) in enumerate(rows):
        tx = templates[(s(d['spectype']), s(d['subtype']).strip('-'))]
        rest = tx.flux.T.dot(d['coeff'][:tx.nbasis]).T
        for b in sp.bands:
            cam[b][j] = sub.R[b][j].dot(resample_flux(sp.wave[b], tx.wave * (1 + d['z']), rest))
    return merge_camera_models(sub, cam, wave), idx

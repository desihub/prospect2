#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.encode
================

Compact "q2" encoding of the flux, noise and model arrays for the page payload
(decoded in the browser by decodeQ2() in the viewer template).
"""

import gzip, base64
import numpy as np


def _b64gz(a):
    #- mtime=0: no timestamp in the gzip header, so the same inputs give a byte-identical page
    return base64.b64encode(gzip.compress(np.ascontiguousarray(a).tobytes(), compresslevel=9, mtime=0)).decode('ascii')


def _esc8(d, masked=None):
    """
    Pack integers as int8 with escapes, for gzip (which compresses bytes far better than int16/32):
      -128 = escape: the value is the next entry of a separate int32 stream
      -127 = masked (only if `masked` is given; then unmasked -127 values are escaped too)
    Returns [int8 stream, int32 escape stream], each gzip+base64.  Lossless.
    """
    d = np.asarray(d, dtype=np.int64)
    lo = -126 if masked is not None else -127
    esc = (d < lo) | (d > 127)
    if masked is not None:
        esc &= ~masked
    a = np.where(esc, -128, d)
    if masked is not None:
        a = np.where(masked, -127, a)
    big = d[esc]
    assert np.all(np.abs(big) < 2**31), 'value too large for int32 escape'
    return [_b64gz(a.astype('i1')), _b64gz(big.astype('<i4'))]


def encode(flux, noise, model, model_precision=100):
    """
    Visually-lossless compact encoding of [nspec, nwave] arrays ("q2"):
      noise : round(100*ln(noise)), delta-coded along wavelength (1% steps)
      flux  : round(10*flux/noise), i.e. quantized at 0.1 sigma against the decoded noise; masked flagged
      model : round(model/scale_i) delta-coded, scale_i = median(noise_i)/model_precision;
              for a spectrum with no unmasked pixels, scale_i = span(model_i)/(10*model_precision)
              where span = max(model,0)-min(model,0).  The model is noiseless, so coarse steps
              show as staircases: 1000 is visually lossless even zoomed+smoothed; 100 can show
              faint few-pixel steps in zoomed, heavily smoothed flat regions.
    Each is packed by _esc8 (int8 + int32 escapes; same values as the old int16/int32 "q1"
    encoding, ~20% smaller after gzip) and decoded in JS by decodeQ2().
    """
    bad = ~np.isfinite(flux) | ~(noise > 0)
    ln = np.round(100 * np.log(np.where(bad, 1.0, noise))).astype(np.int64)
    #- carry previous valid noise into masked pixels so deltas stay ~0
    for i in range(ln.shape[0]):
        idx = np.where(~bad[i], np.arange(ln.shape[1]), 0)
        np.maximum.accumulate(idx, out=idx)
        ln[i] = ln[i][idx]
    dln = np.diff(ln, axis=1, prepend=0)
    noise_q = np.exp(ln / 100)    #- what the decoder will reconstruct
    q = np.round(np.where(bad, 0, flux / noise_q) * 10).astype(np.int64)
    out = dict(enc='q2', noise=_esc8(dln), flux=_esc8(q, masked=bad), model=None, mscale=None)
    if model is not None:
        out['model'], out['mscale'] = encode_model(model, noise, model_precision)
    return out


def encode_model(model, noise, model_precision=100):
    """
    Model rows [nrow, nwave] in steps of scale_i = median(noise_i)/model_precision (noise_i = the
    data noise of the spectrum the row belongs to, NaN where masked), delta-coded and packed by
    _esc8.  See encode().  Returns (packed, per-row scales).
    """
    bad = ~(noise > 0)
    with np.errstate(all='ignore'):
        scale = np.array([np.nanmedian(np.where(b, np.nan, n)) / model_precision if (~b).any() else np.nan
                          for b, n in zip(bad, noise)])
    span = np.ptp(np.concatenate([np.nan_to_num(model), np.zeros((len(model), 1))], axis=1), axis=1)
    fallback = span / (10 * model_precision)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale,
                     np.where(fallback > 0, fallback, 1.0))    #- 1.0: model identically 0
    #- round the stored scale first and quantize with it, so the decoder's error stays <= scale/2
    scale = np.array([float(f'{s:.6g}') for s in scale])
    mq = np.round(np.nan_to_num(model) / scale[:, None]).astype(np.int64)
    return _esc8(np.diff(mq, axis=1, prepend=0)), [float(s) for s in scale]

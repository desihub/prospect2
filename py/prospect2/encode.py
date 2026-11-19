#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.encode
================

Compact "q2" encoding of the flux, noise and model arrays for the page payload
(decoded in the browser by decodeQ2() in the viewer template).

The steps are kept separate so that pages can be sliced and stacked losslessly
at the integer level (prospect2.pagedata):

1. quantize: floats -> integers (quantize, quantize_model)
2. esc8: integers -> int8 matrix [nrow, nwave] + int32 escape values in row order
3. pack: gzip + base64 of both (pack); unpack reverses it
"""

import gzip, base64
import numpy as np


#- gzip level 6, not 9: 7x faster for +1% page size (doc/implementation.md "Payload encoding")
GZIP_LEVEL = 6


def _b64gz(a):
    #- mtime=0: no timestamp in the gzip header, so the same inputs give a byte-identical page
    return base64.b64encode(gzip.compress(np.ascontiguousarray(a).tobytes(), compresslevel=GZIP_LEVEL,
                                          mtime=0)).decode('ascii')


def _unb64gz(b64, dtype):
    return np.frombuffer(gzip.decompress(base64.b64decode(b64)), dtype=dtype)


def esc8(d, masked=None):
    """
    Split integers `d` into an int8 array (same shape) and an int32 escape array:
      -128 = escape: the value is the next entry of the escape array (row-major order)
      -127 = masked (only if `masked` is given; then unmasked -127 values are escaped too)
    gzip compresses the int8 bytes far better than int16/32.  Lossless.
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
    return a.astype('i1'), big.astype('<i4')


def pack(i8, esc):
    """[int8 stream, int32 escape stream], each gzip+base64 (the payload's "data")"""
    return [_b64gz(np.asarray(i8, dtype='i1')), _b64gz(np.asarray(esc, dtype='<i4'))]


def unpack(data, nwave):
    """Inverse of pack: (int8 [nrow, nwave], int32 escapes)"""
    return _unb64gz(data[0], 'i1').reshape(-1, nwave), _unb64gz(data[1], '<i4')


def _esc8(d, masked=None):
    """pack(*esc8(d, masked)): [int8 stream, int32 escape stream], each gzip+base64"""
    return pack(*esc8(d, masked))


def unesc8(i8, esc):
    """Integers from an esc8 (int8, escapes) pair: int64, same shape as i8 (-127 kept as is)"""
    out = np.asarray(i8, dtype=np.int64).copy()
    out[np.asarray(i8) == -128] = esc
    return out


def quantize(flux, noise, model=None, model_precision=100):
    """
    Visually-lossless quantization of [nspec, nwave] arrays ("q2").  Returns
    {name: dict(i8=, esc=[, scale=])} for noise, flux and (if given) model:
      noise : round(100*ln(noise)), delta-coded along wavelength (1% steps)
      flux  : round(10*flux/noise), i.e. quantized at 0.1 sigma against the decoded noise; masked flagged
      model : see quantize_model
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
    out = {}
    out['noise'] = dict(zip(('i8', 'esc'), esc8(dln)))
    out['flux'] = dict(zip(('i8', 'esc'), esc8(q, masked=bad)))
    if model is not None:
        out['model'] = quantize_model(model, noise, model_precision)
    return out


def quantize_model(model, noise, model_precision=100):
    """
    Model rows [nrow, nwave] in steps of scale_i = median(noise_i)/model_precision (noise_i = the
    data noise of the spectrum the row belongs to, NaN where masked), delta-coded.  For a row
    whose spectrum has no unmasked pixels, scale_i = span(model_i)/(10*model_precision), where
    span = max(model,0)-min(model,0).  The model is noiseless, so coarse steps show as
    staircases: 1000 is visually lossless even zoomed+smoothed; 100 can show faint few-pixel
    steps in zoomed, heavily smoothed flat regions.  Returns dict(i8=, esc=, scale=).
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
    i8, esc = esc8(np.diff(mq, axis=1, prepend=0))
    return dict(i8=i8, esc=esc, scale=scale)


#- payload metadata of each q2 array: its row space and encoding
ARRAY_INFO = {
    'noise': dict(rows='spec', enc='q2-dlnoise'),
    'flux': dict(rows='spec', enc='q2-flux'),
    'model': dict(rows='spec', enc='q2-model'),
    'amodel': dict(rows='alt', enc='q2-model'),
}


def decode_noise_flux(noise, flux):
    """
    Floats from the quantized noise and flux (dicts with i8, esc): (flux,
    noise) [nspec, nwave], NaN where masked.  Mirrors decodeQ2() in the viewer.
    """
    ln = np.cumsum(unesc8(noise['i8'], noise['esc']), axis=1)
    sigma = np.exp(ln / 100)
    masked = (np.asarray(flux['i8']) == -127)
    f = unesc8(flux['i8'], flux['esc']) * sigma / 10
    f[masked] = np.nan
    sigma[masked] = np.nan
    return f, sigma


def decode_model(m):
    """Floats [nrow, nwave] from a quantized model (dict with i8, esc, scale)"""
    return np.cumsum(unesc8(m['i8'], m['esc']), axis=1) * np.asarray(m['scale'])[:, None]

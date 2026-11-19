#
# See top-level LICENSE.rst file for Copyright information

"""Tier 1: "q2" payload encoding round trips (decoded by the Python mirror of decodeQ2)"""

import numpy as np

from prospect2.encode import quantize, pack, unpack, _esc8, decode_noise_flux, decode_model
from prospect2.test.util import decode_q2, _unesc8


def test_esc8_roundtrip():
    d = np.array([0, 1, -1, 127, -126, -127, -128, 128, -129, 10**6, -10**6])
    out, _ = _unesc8(_esc8(d))
    assert np.array_equal(out, d)
    #- with a mask: masked entries become -127, and genuine -127 values are escaped
    masked = np.zeros(len(d), dtype=bool)
    masked[[0, 3]] = True
    out, m = _unesc8(_esc8(d, masked=masked), masked_flag=True)
    assert np.array_equal(m, masked)
    assert np.array_equal(out[~masked], d[~masked])


def _payload(flux, noise, model, precision=100):
    """the payload's q2 part (nspec, nwave, arrays) for these arrays"""
    nspec, nwave = flux.shape
    arrays = {}
    for name, q in quantize(flux, noise, model, model_precision=precision).items():
        arrays[name] = dict(data=pack(q['i8'], q['esc']))
        if 'scale' in q:
            arrays[name]['scale'] = [float(x) for x in q['scale']]
    return dict(nspec=nspec, nwave=nwave, arrays=arrays)


def test_encode_roundtrip():
    rng = np.random.default_rng(2)
    nspec, nwave = 5, 3000
    noise = np.exp(rng.normal(0, 1, (nspec, nwave)))
    flux = rng.normal(0, 1, (nspec, nwave)) * noise * 3
    flux[0, 100] = 1e4 * noise[0, 100]             #- far beyond int8: escaped
    flux[1, 200:300] = np.nan                      #- masked
    noise[1, 200:300] = np.nan
    flux[2, :] = np.nan                            #- fully masked spectrum
    noise[2, :] = np.nan
    model = np.cumsum(rng.normal(0, 1, (nspec, nwave)), axis=1)
    p = _payload(flux, noise, model)
    dflux, dnoise, dmodel, damodel = decode_q2(p)

    bad = ~np.isfinite(flux)
    assert np.array_equal(np.isnan(dflux), bad), 'mask not preserved'
    good = ~bad
    #- noise in 1% log steps; flux quantized at 0.1 sigma of the decoded noise
    assert np.all(np.abs(np.log(dnoise[good] / noise[good])) <= 0.005 + 1e-12)
    assert np.all(np.abs(dflux[good] - flux[good]) <= 0.05 * dnoise[good] * (1 + 1e-9))
    #- model: steps of median(noise)/precision, error <= half a step; fallback scale if no data
    scale = np.array(p['arrays']['model']['scale'])
    assert np.all(np.abs(dmodel - model) <= 0.5 * scale[:, None] * (1 + 1e-9) + 1e-12 * np.abs(model))
    assert np.isclose(scale[0], np.nanmedian(noise[0]) / 100, rtol=1e-5)
    span = max(model[2].max(), 0) - min(model[2].min(), 0)
    assert np.isclose(scale[2], span / 1000, rtol=1e-5)
    assert damodel is None

    #- the package's own decoder (PageData.flux etc.) agrees with the JS mirror
    q = quantize(flux, noise, model)
    pflux, pnoise = decode_noise_flux(q['noise'], q['flux'])
    assert np.array_equal(pflux, dflux, equal_nan=True) and np.array_equal(pnoise, dnoise, equal_nan=True)
    assert np.array_equal(decode_model(q['model']), dmodel)
    #- pack / unpack are inverses
    i8, esc = unpack(pack(q['flux']['i8'], q['flux']['esc']), nwave)
    assert np.array_equal(i8, q['flux']['i8']) and np.array_equal(esc, q['flux']['esc'])


def test_encode_no_model():
    flux = np.ones((2, 10))
    noise = np.ones((2, 10))
    p = _payload(flux, noise, None)
    assert 'model' not in p['arrays']
    dflux, dnoise, dmodel, _ = decode_q2(p)
    assert np.allclose(dflux, 1) and dmodel is None

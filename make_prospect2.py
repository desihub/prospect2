#!/usr/bin/env python

"""
Prototype "prospect2" spectral viewer generator.

Writes a single self-contained HTML file: a small custom JS app (uPlot for the
main plot, plain canvas for thumbnails) plus a compact gzip+base64 float16
data payload.  No Bokeh, no server, works offline.

Usage:
    python alt/make_prospect2.py -i data/coadd-4-1496-thru20240214.fits \
        -o alt/prospect2_example.html -f 2000:2020
"""

import os, sys, json, gzip, base64, argparse, warnings
import numpy as np
from astropy.table import Table
from astropy.io import fits

from desispec.io import read_spectra
from desispec.io.util import checkgzip
from desispec.coaddition import coadd_cameras

HERE = os.path.dirname(os.path.abspath(__file__))


def _b64gz(a):
    return base64.b64encode(gzip.compress(np.ascontiguousarray(a).tobytes(), compresslevel=9, mtime=0)).decode('ascii')


def encode(flux, noise, model, model_precision=100):
    """
    Visually-lossless compact encoding of [nspec, nwave] arrays:
      noise : round(100*ln(noise)) int16, delta-coded along wavelength (1% steps)
      flux  : round(10*flux/noise) int16, i.e. quantized at 0.1 sigma; -32768 = masked
      model : round(model/scale_i) int32 delta-coded, scale_i = median(noise_i)/model_precision;
              for a spectrum with no unmasked pixels, scale_i = span(model_i)/(10*model_precision)
              where span = max(model,0)-min(model,0).  The model is noiseless, so coarse steps
              show as staircases: 1000 is visually lossless even zoomed+smoothed; 100 can show
              faint few-pixel steps in zoomed, heavily smoothed flat regions.
    Each is gzip+base64.  Decoded in JS by decodeQ().
    """
    bad = ~np.isfinite(flux) | ~(noise > 0)
    ln = np.round(100 * np.log(np.where(bad, 1.0, noise))).astype(np.int32)
    #- carry previous valid noise into masked pixels so deltas stay ~0
    for i in range(ln.shape[0]):
        idx = np.where(~bad[i], np.arange(ln.shape[1]), 0)
        np.maximum.accumulate(idx, out=idx)
        ln[i] = ln[i][idx]
    dln = np.diff(ln, axis=1, prepend=0).astype('<i2')
    noise_q = np.exp(ln / 100)    #- what the decoder will reconstruct
    q = np.round(np.where(bad, 0, flux / noise_q) * 10)
    qtype = '<i2' if np.all(np.abs(q) < 32767) else '<i4'    #- int32 only if S/N > 3276/pix
    q = q.astype(qtype)
    q[bad] = np.iinfo(qtype).min
    out = dict(enc='q1', noise=_b64gz(dln), flux=_b64gz(q), fluxbits=8*np.dtype(qtype).itemsize,
               model=None, mscale=None)
    if model is not None:
        with np.errstate(all='ignore'):
            scale = np.array([np.nanmedian(np.where(b, np.nan, n)) / model_precision if (~b).any() else np.nan
                              for b, n in zip(bad, noise)])
        span = np.ptp(np.concatenate([np.nan_to_num(model), np.zeros((len(model), 1))], axis=1), axis=1)
        fallback = span / (10 * model_precision)
        scale = np.where(np.isfinite(scale) & (scale > 0), scale,
                         np.where(fallback > 0, fallback, 1.0))    #- 1.0: model identically 0
        mq = np.round(np.nan_to_num(model) / scale[:, None]).astype(np.int64)
        dmq = np.diff(mq, axis=1, prepend=0)
        assert np.abs(dmq).max() < 2**31, 'model too large for int32 delta encoding'
        out['model'] = _b64gz(dmq.astype('<i4'))
        out['mscale'] = [float(f'{s:.6g}') for s in scale]
    return out


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


def find_file(filename):
    """filename or filename(.gz) if it exists, else None"""
    try:
        return checkgzip(filename)
    except (FileNotFoundError, IOError):
        return None


def read_rrmodel(rrmodel_file, sp):
    """
    Read per-camera redrock models {band: [nspec, nwave_band]} from an rrmodel
    file (REDSHIFTS + {B,R,Z}_WAVELENGTH + {B,R,Z}_MODEL HDUs), matched to the
    rows of sp by TARGETID.  Returns None, with a warning, if the file doesn't
    cover every spectrum and camera on the same wavelength grids.
    """
    with fits.open(rrmodel_file, memmap=True) as hx:
        names = {h.name for h in hx}
        tids = np.asarray(hx['REDSHIFTS'].data['TARGETID'])
        row = {t: i for i, t in enumerate(tids)}
        want = np.asarray(sp.fibermap['TARGETID'])
        missing = [t for t in want if t not in row]
        if missing:
            print(f'WARNING: {len(missing)} TARGETIDs not in {rrmodel_file}; computing models instead')
            return None
        idx = np.array([row[t] for t in want])
        out = {}
        for b in sp.bands:
            B = b.upper()
            if f'{B}_MODEL' not in names or f'{B}_WAVELENGTH' not in names:
                print(f'WARNING: no {B}_MODEL in {rrmodel_file}; computing models instead')
                return None
            w = hx[f'{B}_WAVELENGTH'].data
            if len(w) != len(sp.wave[b]) or not np.allclose(w, sp.wave[b]):
                print(f'WARNING: {B}_WAVELENGTH in {rrmodel_file} differs from spectra; computing models instead')
                return None
            out[b] = hx[f'{B}_MODEL'].data[idx].astype(np.float64)
    return out


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


def parse_fiberlist(spec):
    """
    "10,20,30-40,50:60" -> [10, 20, 30..40, 50..59]: "a-b" is inclusive,
    "a:b" is python-style (excludes b).  Same syntax as the viewer's filter box.
    """
    fibers = []
    for part in spec.replace(' ', '').split(','):
        if not part:
            continue
        if ':' in part:
            a, b = part.split(':'); fibers.extend(range(int(a), int(b)))
        elif '-' in part[1:]:
            a, b = part.split('-'); fibers.extend(range(int(a), int(b) + 1))
        else:
            fibers.append(int(part))
    return np.array(fibers, dtype=int)


def air_to_vac(air):
    s2 = (1e4 / air) ** 2
    n = 1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - s2) \
        + 1.599740894897e-4 / (3.892568793293e1 - s2)
    return air * n


def read_lines():
    lines = []
    for kind in ('emission', 'absorption'):
        fn = os.path.join(HERE, 'data', f'{kind}_lines.txt')
        for row in open(fn, encoding='utf-8'):
            if row.startswith('#') or not row.strip():
                continue
            name, longname, lam, vac, major = [x.strip() for x in row.split(',')]
            lam = float(lam)
            if vac != 'True':
                lam = air_to_vac(lam)
            lines.append(dict(name=name, longname=longname, wave=round(lam, 2),
                              kind=kind[:3], major=(major == 'True')))
    return lines


def mag(flux):
    flux = float(flux)
    return round(22.5 - 2.5 * np.log10(flux), 2) if flux > 0 else None


def target_classes(row):
    try:
        from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
    except ImportError:
        return [], []
    names = []
    for col, mask in [('DESI_TARGET', desi_mask), ('BGS_TARGET', bgs_mask),
                      ('MWS_TARGET', mws_mask), ('SCND_TARGET', scnd_mask)]:
        if col in row.colnames and row[col] != 0:
            names += [n for n in mask.names(int(row[col]))
                      if n not in ('BGS_ANY', 'MWS_ANY', 'SCND_ANY') and not n.endswith(('_SOUTH', '_NORTH'))]
    #- primary classes for compact badges
    primary = [n for n in names if n in ('ELG', 'LRG', 'QSO', 'BGS_BRIGHT', 'BGS_FAINT', 'MWS_BROAD',
                                         'MWS_MAIN_BLUE', 'MWS_MAIN_RED', 'STD_FAINT', 'STD_BRIGHT',
                                         'STD_WD', 'MWS_WD', 'ELG_LOP', 'ELG_VLO', 'ELG_HIP')]
    return names, primary


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--infile', required=True, help='input coadd/spectra file')
    p.add_argument('-o', '--outfile', required=True, help='output html file')
    p.add_argument('-f', '--fibers', help='fibers, e.g. 10,20,30-40,50:60 (a-b inclusive, a:b excludes b)')
    p.add_argument('-t', '--targetids', help='comma-separated TARGETIDs (combined with -f as a union)')
    p.add_argument('--title', default=None)
    p.add_argument('--no-rrmodel', action='store_true',
                   help='ignore any rrmodel*.fits file; compute models from redrock templates + COEFF')
    p.add_argument('--model-precision', type=float, default=100,
                   help='store the model in steps of (median data noise)/N; default %(default)s. '
                        '1000 is visually lossless even zoomed+smoothed (~+1.6 MB per 500 spectra)')
    args = p.parse_args()

    sp = read_spectra(args.infile)

    basedir, basefile = os.path.split(args.infile)
    rrfile = basefile.replace('coadd', 'redrock').replace('spectra', 'redrock')
    rrfile = find_file(os.path.join(basedir, rrfile))
    zcat = None
    rrmodel_file = None
    if rrfile:
        zcat = Table.read(rrfile, 'REDSHIFTS')
        zcat.meta.update(fits.getheader(rrfile, 0))
        rrdir, rrbase = os.path.split(rrfile)
        if not args.no_rrmodel and rrbase.startswith('redrock'):
            rrmodel_file = find_file(os.path.join(rrdir, 'rrmodel' + rrbase[len('redrock'):]))

    #- -f and -t select the union of the listed fibers and TARGETIDs (all spectra if neither)
    keep = np.zeros(sp.num_spectra(), dtype=bool) if (args.fibers or args.targetids) else np.ones(sp.num_spectra(), dtype=bool)
    if args.fibers:
        keep |= np.isin(sp.fibermap['FIBER'], parse_fiberlist(args.fibers))
    if args.targetids:
        keep |= np.isin(sp.fibermap['TARGETID'], [int(t) for t in args.targetids.replace(' ', '').split(',') if t])
    sp = sp[keep]
    if zcat is not None:
        zcat = zcat[keep]
    if sp.num_spectra() == 0:
        sys.exit('ERROR: no spectra selected')

    #- Per-camera median S/N (before merging cameras)
    snr = {}
    for b in sp.bands:
        with np.errstate(invalid='ignore'):
            s = sp.flux[b] * np.sqrt(sp.ivar[b])
        with warnings.catch_warnings():         #- fully masked fibers -> NaN (shown as "—")
            warnings.simplefilter('ignore', RuntimeWarning)
            snr[b] = np.nanmedian(np.where(sp.ivar[b] > 0, s, np.nan), axis=1)

    #- Camera overlap regions [lo, hi, label], from the per-camera grids
    bands = sorted(sp.bands, key=lambda b: sp.wave[b][0])
    overlaps = [[float(sp.wave[b2][0]), float(sp.wave[b1][-1]), f'{b1}/{b2} overlap']
                for b1, b2 in zip(bands[:-1], bands[1:]) if sp.wave[b2][0] < sp.wave[b1][-1]]

    #- Redrock model on per-camera grids, then merged
    model = None
    model_source = None
    if zcat is not None:
        cam_models = read_rrmodel(rrmodel_file, sp) if rrmodel_file else None
        if cam_models is not None:
            model_source = os.path.basename(rrmodel_file)
            print(f'INFO: redrock models read from {rrmodel_file}')
        else:
            cam_models = camera_models(sp, zcat)
            model_source = 'redrock templates + COEFF'
            print('INFO: redrock models computed from templates + COEFF')

    #- Merge cameras onto a single linear grid
    sp_all = coadd_cameras(sp)
    band = sp_all.bands[0]
    wave = sp_all.wave[band]
    flux = sp_all.flux[band].copy()
    ivar = sp_all.ivar[band]
    bad = (ivar <= 0)
    if sp_all.mask is not None:
        bad |= (sp_all.mask[band] != 0)
    flux[bad] = np.nan
    with np.errstate(divide='ignore'):
        noise = np.where(bad, np.nan, 1 / np.sqrt(ivar))
    if zcat is not None:
        model = merge_camera_models(sp, cam_models, wave)

    dw = np.diff(wave)
    assert np.allclose(dw, dw[0]), 'expected a linear wavelength grid'

    #- Per-target metadata
    from redrock.zwarning import ZWarningMask
    zw_names = [n for n in dir(ZWarningMask) if n.isupper()]
    try:
        from desispec.maskbits import fibermask
    except ImportError:
        fibermask = None

    fm = sp.fibermap
    targets = []
    for i in range(sp.num_spectra()):
        r = fm[i]
        allclass, primary = target_classes(r)
        t = dict(
            i=i, targetid=str(int(r['TARGETID'])), fiber=int(r['FIBER']) if 'FIBER' in fm.colnames else None,
            tileid=int(r['TILEID']) if 'TILEID' in fm.colnames else None,
            ra=round(float(r['TARGET_RA']), 6), dec=round(float(r['TARGET_DEC']), 6),
            objtype=str(r['OBJTYPE']), morph=str(r['MORPHTYPE']).strip(),
            classes=allclass, primary=primary,
            mag={k: mag(r[f'FLUX_{k.upper()}']) for k in ('g', 'r', 'z', 'w1', 'w2')},
            fibermag_r=mag(r['FIBERFLUX_R']),
            ebv=round(float(r['EBV']), 3),
            exptime=round(float(r['COADD_EXPTIME']), 1), numexp=int(r['COADD_NUMEXP']),
            fiberstatus=int(r['COADD_FIBERSTATUS']),
            fiberstatus_names=(fibermask.names(int(r['COADD_FIBERSTATUS'])) if fibermask else []),
            snr={b: (round(float(snr[b][i]), 2) if np.isfinite(snr[b][i]) else None) for b in sp.bands},
            #- cameras with no unmasked pixels for this spectrum: [lo, hi, band]
            masked_cams=[[float(sp.wave[b][0]), float(sp.wave[b][-1]), b] for b in bands
                         if not np.any((sp.ivar[b][i] > 0) & ((sp.mask[b][i] == 0) if sp.mask is not None else True))],
        )
        if zcat is not None:
            zb = zcat[i]
            zwarn = int(zb['ZWARN'])
            t.update(z=float(zb['Z']), zerr=float(zb['ZERR']), zwarn=zwarn,
                     zwarn_names=[n for n in zw_names if zwarn & getattr(ZWarningMask, n)],
                     spectype=str(zb['SPECTYPE']).strip(), subtype=str(zb['SUBTYPE']).strip().strip('-'),
                     deltachi2=round(float(zb['DELTACHI2']), 2), chi2=round(float(zb['CHI2']), 1),
                     npixels=int(zb['NPIXELS']))
        targets.append(t)

    payload = dict(
        title=args.title or f'{basefile}',
        source=basefile,
        redrock=os.path.basename(rrfile) if zcat is not None else None,
        model_source=model_source,
        wave0=float(wave[0]), dwave=float(dw[0]), nwave=len(wave), nspec=sp.num_spectra(),
        flux_unit='1e-17 erg/s/cm²/Å',
        targets=targets,
        lines=read_lines(),
        overlaps=overlaps,
        **encode(flux, noise, model, model_precision=args.model_precision),
    )

    template = open(os.path.join(HERE, 'prospect2_template.html'), encoding='utf-8').read()
    uplot_js = open(os.path.join(HERE, 'vendor', 'uPlot.iife.min.js'), encoding='utf-8').read()
    uplot_css = open(os.path.join(HERE, 'vendor', 'uPlot.min.css'), encoding='utf-8').read()
    html = (template.replace('/*__UPLOT_CSS__*/', uplot_css)
                    .replace('/*__UPLOT_JS__*/', uplot_js)
                    .replace('"__DATA__"', json.dumps(payload, separators=(',', ':')))
                    .replace('__TITLE__', payload['title']))
    with open(args.outfile, 'w', encoding='utf-8') as fx:
        fx.write(html)
    print(f'Wrote {args.outfile}: {os.path.getsize(args.outfile)/1e6:.2f} MB, {sp.num_spectra()} spectra')


if __name__ == '__main__':
    main()

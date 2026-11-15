#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.scripts.plot_prospect2
================================

Command-line interface: read a coadd (or spectra) file and its matching redrock
files, then write one offline viewer page with prospect2.make_html().

Usage (with bin/ on $PATH and py/ on $PYTHONPATH):
    plot_prospect2 -i coadd-4-1496-thru20240214.fits    # -> ./prospect-4-1496-thru20240214.html
    plot_prospect2 -i coadd-4-1496-thru20240214.fits -o outdir/
    plot_prospect2 -i coadd-4-1496-thru20240214.fits \\
        -o example-20.html -f 2000:2020
"""

import os, sys, argparse

from desiutil.log import get_logger
from desispec.io import read_spectra

from prospect2.page import make_html, write_html
from prospect2.io import find_redrock_files, read_redshifts, read_rrmodel, read_rrdetails
from prospect2.util import select_spectra, resolve_outfile


def parse(options=None):
    """
    Parse command-line `options` (list of str; default sys.argv[1:]).  Resolves
    args.outfile to the output file path (see resolve_outfile), exiting with an
    error message if no output name can be derived.
    """
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--infile', required=True, help='input coadd/spectra file')
    p.add_argument('-o', '--outfile', default=None,
                   help='output html file, or directory (existing, or ending in "/"; created if '
                        'needed) for the default name: coadd-X.fits -> prospect-X.html, in the '
                        'current directory if -o is not given.  The default name is only used '
                        'without -f/-t (pages with every input spectrum)')
    p.add_argument('-f', '--fibers', help='fibers, e.g. 10,20,30-40,50:60 (a-b inclusive, a:b excludes b)')
    p.add_argument('-t', '--targetids', help='comma-separated TARGETIDs (combined with -f as a union)')
    p.add_argument('--title', default=None)
    p.add_argument('--no-rrmodel', action='store_true',
                   help='ignore any rrmodel*.fits file; compute models from redrock templates + COEFF')
    p.add_argument('--nfits', type=int, default=5,
                   help='number of best redrock fits to keep from a matching rrdetails*.h5 file (0 = ignore it)')
    p.add_argument('--alt-models', type=int, default=2,
                   help='store real redrock models for up to this many alternate fits per target '
                        '(fits 2..N+1, needs rrdetails; 0 = none, the model is then hidden for alternates)')
    p.add_argument('--alt-max-dchi2', type=float, default=100,
                   help='only store a model for alternate fits with chi2 - best chi2 below this')
    p.add_argument('--alt-model-precision', type=float, default=30,
                   help='as --model-precision, for alternate-fit models')
    p.add_argument('--model-precision', type=float, default=100,
                   help='store the model in steps of (median data noise)/N; default %(default)s. '
                        '1000 is visually lossless even zoomed+smoothed (~+1.2 MB per 500 spectra)')
    args = p.parse_args(options)

    #- resolve the output before the slow read, so naming errors are immediate
    try:
        args.outfile = resolve_outfile(args.outfile, args.infile, subset=bool(args.fibers or args.targetids))
    except ValueError as err:
        p.error(str(err))
    return args


def main(args):
    log = get_logger()
    outdir = os.path.dirname(args.outfile)
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir)
        log.info(f'Created directory {outdir}')

    sp = read_spectra(args.infile)
    rr = find_redrock_files(args.infile)
    zcat = read_redshifts(rr['redrock']) if rr['redrock'] else None

    keep = select_spectra(sp.fibermap, args.fibers, args.targetids)
    if not keep.any():
        log.error('no spectra selected')
        return 1
    sp = sp[keep]
    if zcat is not None:
        zcat = zcat[keep]

    models = fits = None
    use_rrmodel = zcat is not None and rr['rrmodel'] and not args.no_rrmodel
    if use_rrmodel:
        models = read_rrmodel(rr['rrmodel'], sp.fibermap['TARGETID'], sp.wave)
        log.info(f'redrock models read from {rr["rrmodel"]}')
    if zcat is not None and rr['rrdetails'] and args.nfits > 0:
        fits = read_rrdetails(rr['rrdetails'], zcat, nfits=args.nfits)
        log.info(f'top {args.nfits} redrock fits read from {rr["rrdetails"]}')

    base = os.path.basename
    provenance = dict(source=base(args.infile),
                      redrock=base(rr['redrock']) if rr['redrock'] else None,
                      model_source=base(rr['rrmodel']) if use_rrmodel else None,
                      rrdetails=base(rr['rrdetails']) if fits is not None else None)
    html = make_html(sp, zcat=zcat, models=models, fits=fits, title=args.title, provenance=provenance,
                     alt_models=args.alt_models, alt_max_dchi2=args.alt_max_dchi2,
                     model_precision=args.model_precision, alt_model_precision=args.alt_model_precision)
    write_html(html, args.outfile)
    log.info(f'{sp.num_spectra()} spectra')
    return 0


if __name__ == '__main__':
    sys.exit(main(parse()))

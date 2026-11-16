#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.util
==============

Spectrum selection and output naming.
"""

import os
import numpy as np

from desiutil.log import get_logger


def parse_fiberlist(spec):
    """
    "10,20,30-40,50:60" -> [10, 20, 30..40, 50..59]: "a-b" is inclusive,
    "a:b" is python-style (excludes b).  Same syntax as the viewer's filter box
    (parseIdList() in the template; the browser test checks they agree), except
    that the filter box reads numbers >= 5000 or < 0 as TARGETIDs.
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


def fibermap_value(fibermap, i, *names, default=None):
    """
    Row `i` of the first of the columns `names` that `fibermap` has, else
    `default` (also for a masked value).  Only TARGETID is required to make a
    page; every other fibermap column it reads is optional, via this function.
    """
    for name in names:
        if name in fibermap.colnames:
            value = fibermap[name][i]
            return default if np.ma.is_masked(value) else value
    return default


def select_spectra(fibermap, fibers=None, targetids=None):
    """
    Boolean mask of the fibermap rows to keep: the union of `fibers` (a list
    string as in parse_fiberlist) and `targetids` (comma-separated string), or
    every row if neither is given.
    """
    if not (fibers or targetids):
        return np.ones(len(fibermap), dtype=bool)
    keep = np.zeros(len(fibermap), dtype=bool)
    if fibers:
        if 'FIBER' not in fibermap.colnames:
            raise ValueError('these spectra have no FIBER column (e.g. a healpix coadd); select with -t TARGETIDs')
        keep |=np.isin(fibermap['FIBER'], parse_fiberlist(fibers))
    if targetids:
        keep |= np.isin(fibermap['TARGETID'], [int(t) for t in targetids.replace(' ', '').split(',') if t])
    return keep


def default_outfile(infile):
    """
    Default output basename for a page holding every spectrum in `infile`:
    coadd-X.fits[.gz] or spectra-X.fits[.gz] -> prospect-X.html.  None if the
    input name doesn't follow that pattern.
    """
    base = os.path.basename(infile)
    for ext in ('.fits.gz', '.fits'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    else:
        return None
    for prefix in ('coadd-', 'spectra-'):
        if base.startswith(prefix) and len(base) > len(prefix):
            return 'prospect-' + base[len(prefix):] + '.html'
    return None


def resolve_outfile(outarg, infile, subset):
    """
    Output file path from the -o argument (None if not given).  cp-style: an
    existing directory or a path ending in '/' is a directory, which gets the
    default name; anything else is a file name.  Raises ValueError with a message
    explaining why a file name is needed if the default name can't be used.
    Warns if an explicit file name doesn't end in .html/.htm.
    """
    if outarg is not None and not (outarg.endswith(os.sep) or os.path.isdir(outarg)):
        if not outarg.lower().endswith(('.html', '.htm')):
            get_logger().warning(f"output {outarg} doesn't end in .html; writing it as a file "
                                 f"(use {outarg}/ for a directory)")
        return outarg

    where = 'no -o was given' if outarg is None else f'-o {outarg} is a directory'
    name = default_outfile(infile)
    if subset:
        dflt = f' {name}' if name else ''
        raise ValueError(f'{where}, but -f/-t select a subset of the spectra in '
                         f'{os.path.basename(infile)}, so the default name{dflt} (reserved for pages '
                         f'with every spectrum in the input) does not apply; give -o FILE.html')
    if name is None:
        raise ValueError(f"{where}, but can't derive an output name from {os.path.basename(infile)}: "
                         f"default names map coadd-*.fits or spectra-*.fits to prospect-*.html; "
                         f"give -o FILE.html")
    return os.path.join(outarg, name) if outarg else name

#
# See top-level LICENSE.rst file for Copyright information

"""
prospect2.grouping
==================

How the spectra on a page are grouped (header keyword SPGRP), and which
columns uniquely identify a spectrum.  The viewer keys visual inspection (VI)
on that key and exports its columns.

* tile-based (cumulative, pernight, perexp): TILEID, FIBER, NIGHT, plus EXPID
  for uncoadded (per-exposure) spectra.  A TARGETID can be on several tiles.
* cross-tile (healpix, uniqpix): SURVEY, PROGRAM, TARGETID.  A TARGETID is
  coadded across tiles within a survey and program, never across them.

The values come from each input's headers, which desispec.spectra.stack keeps
only for the first input, so they are read per input, before stacking.
"""

import numpy as np
from desiutil.log import get_logger

TILE_GROUPS = ('cumulative', 'pernight', 'perexp')
HEALPIX_GROUPS = ('healpix', 'uniqpix')

#- per-row values; KEYCOLS are the possible key columns, in display/export order
COLUMNS = ('SPGRP', 'SURVEY', 'PROGRAM', 'TILEID', 'FIBER', 'NIGHT', 'EXPID')
TILE_KEY = ('TILEID', 'FIBER', 'NIGHT', 'EXPID')
HEALPIX_KEY = ('SURVEY', 'PROGRAM', 'TARGETID')


def _meta(d, name):
    """d[name] for a header-like dict `d` (may be None), else None"""
    return None if d is None else d.get(name)


def _str(x):
    return None if x is None else str(x).strip().lower() or None


def _int(x):
    return None if x is None else int(x)


def grouping_columns(sp):
    """
    Per-row grouping values for one input Spectra (not a stack of several
    inputs): {column: list of nspec values} for COLUMNS, None where unknown.

    Fibermap columns win over headers.  PROGRAM is the primary header's (the
    healpix grouping's, e.g. 'dark'), else the fibermap header's FAPRGRM (the
    tile's, e.g. 'dark1b'), else its PROGRAM.  NIGHT for a cumulative coadd is
    its LASTNIGHT (header NIGHT, or SPGRPVAL).
    """
    n = sp.num_spectra()
    fm = sp.fibermap
    fmeta = fm.meta if fm is not None else None
    spgrp = _str(_meta(sp.meta, 'SPGRP'))

    def column(name, *headers, conv):
        if fm is not None and name in fm.colnames:
            return [conv(v) for v in fm[name]]
        for d, key in headers:
            v = _meta(d, key)
            if v is not None:
                return [conv(v)] * n
        return [None] * n

    spgrpval = (sp.meta, 'SPGRPVAL')
    return dict(
        SPGRP=[spgrp] * n,
        SURVEY=column('SURVEY', (sp.meta, 'SURVEY'), (fmeta, 'SURVEY'), conv=_str),
        PROGRAM=column('PROGRAM', (sp.meta, 'PROGRAM'), (fmeta, 'FAPRGRM'), (fmeta, 'PROGRAM'), conv=_str),
        TILEID=column('TILEID', (sp.meta, 'TILEID'), (fmeta, 'TILEID'), conv=_int),
        FIBER=column('FIBER', conv=_int),
        NIGHT=column('NIGHT', (sp.meta, 'NIGHT'), *([spgrpval] if spgrp in ('cumulative', 'pernight') else []), conv=_int),
        EXPID=column('EXPID', (sp.meta, 'EXPID'), *([spgrpval] if spgrp == 'perexp' else []), conv=_int),
    )


def concat_columns(cols):
    """Concatenate grouping_columns() of several inputs, in order"""
    return {c: [v for d in cols for v in d[c]] for c in COLUMNS}


def key_columns(groups, targetids):
    """
    (spgrp, keycols) for the grouping values `groups` (from grouping_columns,
    possibly concatenated) and the matching TARGETIDs.  spgrp is the common
    SPGRP or None; keycols are the columns that make each spectrum unique,
    a subset of TILE_KEY or HEALPIX_KEY, or ['ROW'] (row index on the page) if
    no subset is unique (logged).

    Raises ValueError if the inputs mix SPGRP groupings.
    """
    log = get_logger()
    spgrps = sorted({g for g in groups['SPGRP'] if g is not None})
    if len(spgrps) > 1:
        raise ValueError(f'inputs mix spectral groupings (SPGRP): {", ".join(spgrps)}')
    spgrp = spgrps[0] if spgrps else None
    vals = dict(groups, TARGETID=[int(t) for t in targetids])
    known = lambda c: all(v is not None for v in vals[c])
    if spgrp in HEALPIX_GROUPS or (spgrp not in TILE_GROUPS and not (known('TILEID') and known('FIBER'))):
        keycols = [c for c in HEALPIX_KEY if known(c)]
    else:
        keycols = [c for c in TILE_KEY if known(c)]
        if keycols[:2] != ['TILEID', 'FIBER']:
            keycols = []
    rows = list(zip(*(vals[c] for c in keycols)))
    if not keycols or len(set(rows)) != len(rows):
        log.warning(f'spectra are not unique by {"/".join(keycols) or "TILEID/FIBER or TARGETID"} '
                    f'(SPGRP={spgrp}); VI is keyed by row on this page')
        keycols = ['ROW']
    return spgrp, keycols

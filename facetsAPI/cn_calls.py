"""Copy-number call table keyed on (wgd, tcn, lcn).

``lcn`` is the LESSER copy number (facets ``lcn.em``). ``None`` stands for an
undetermined lcn (NA in facets output). Entries where lcn is impossible for the
given tcn (lcn > tcn/2) are deliberately absent.

Historical note: the table shipped before 0.6.0 was keyed on major-cn-shaped
values (e.g. ``(False, 1, 1)`` for HETLOSS) while the lookup passed lcn, so most
real segments resolved to ``None``. This module carries the corrected table.
"""

import math

CN_CALL_TABLE = {
    # ------------------ No genome doubling (wgd=False) ------------------
    (False, 0, 0):    'HOMDEL',
    (False, 0, None): 'HOMDEL',

    # HETLOSS
    (False, 1, 0):    'HETLOSS',
    (False, 1, None): 'HETLOSS',

    # CNLOH / DIPLOID (and ambiguous)
    (False, 2, 0):    'CNLOH',
    (False, 2, 1):    'DIPLOID',
    (False, 2, None): None,  # DIPLOID or CNLOH ambiguous

    # GAIN / CNLOH & GAIN
    (False, 3, 0):    'CNLOH & GAIN',
    (False, 3, 1):    'GAIN',
    (False, 3, None): 'GAIN',

    (False, 4, 0):    'CNLOH & GAIN',
    (False, 4, 1):    'GAIN',
    (False, 4, 2):    'TETRAPLOID',
    (False, 4, None): 'GAIN',

    # AMP (and variants)
    (False, 5, 0):    'AMP (LOH)',
    (False, 5, 1):    'AMP',   # (mcn=4, lcn=1)
    (False, 5, 2):    'AMP',   # (mcn=3, lcn=2)
    (False, 5, None): 'AMP',

    (False, 6, 0):    'AMP (LOH)',
    (False, 6, 1):    'AMP',   # (mcn=5, lcn=1)
    (False, 6, 2):    'AMP',   # (mcn=4, lcn=2)
    (False, 6, 3):    'AMP (BALANCED)',  # (mcn=3, lcn=3)
    (False, 6, None): 'AMP',

    (False, 7, 1):    'AMP',   # (mcn=6, lcn=1)
    (False, 7, 2):    'AMP',   # (mcn=5, lcn=2)
    (False, 7, 3):    'AMP',   # (mcn=4, lcn=3)

    (False, 8, 2):    'AMP',   # (mcn=6, lcn=2)
    (False, 8, 3):    'AMP',   # (mcn=5, lcn=3)

    (False, 9, 3):    'AMP',   # (mcn=6, lcn=3)

    # ------------------ With genome doubling (wgd=True) ------------------
    (True, 0, 0):     'HOMDEL',
    (True, 0, None):  'HOMDEL',

    (True, 1, 0):     'LOSS BEFORE & AFTER',
    (True, 1, None):  'LOSS BEFORE & AFTER',

    (True, 2, 0):     'LOSS BEFORE',
    (True, 2, 1):     'DOUBLE LOSS AFTER',
    (True, 2, None):  'LOSS BEFORE or DOUBLE LOSS AFTER',

    (True, 3, 0):     'CNLOH BEFORE & LOSS',
    (True, 3, 1):     'LOSS AFTER',
    (True, 3, None):  'LOSS',

    (True, 4, 0):     'CNLOH BEFORE',
    (True, 4, 1):     'CNLOH AFTER',
    (True, 4, 2):     'TETRAPLOID',
    (True, 4, None):  'TETRAPLOID or CNLOH BEFORE',

    (True, 5, 0):     'CNLOH BEFORE & GAIN',
    (True, 5, 1):     'LOSS & GAIN',
    (True, 5, 2):     'GAIN',
    (True, 5, None):  'GAIN',

    (True, 6, 0):     'AMP (LOH)',
    (True, 6, 1):     'AMP',
    (True, 6, 2):     'AMP',
    (True, 6, 3):     'AMP (BALANCED)',
    (True, 6, None):  'AMP',

    (True, 7, 1):     'AMP',
    (True, 7, 2):     'AMP',
    (True, 7, 3):     'AMP',

    (True, 8, 2):     'AMP',
    (True, 8, 3):     'AMP',

    (True, 9, 3):     'AMP',
}


def _norm_cn(value):
    """Normalize a copy-number value: None/NaN/-1 -> None, otherwise int."""
    if value is None:
        return None
    try:
        if isinstance(value, float) and math.isnan(value):
            return None
        iv = int(value)
    except (TypeError, ValueError):
        return None
    if iv < 0:          # facetsAPI stores NA as -1 on FacetsSegment
        return None
    return iv


def get_cn_call(wgd, tcn, lcn):
    """Return the call string for (wgd, tcn, lcn), or None when the table has no entry.

    Accepts floats, NaN and the -1 NA sentinel used by FacetsSegment.
    """
    t = _norm_cn(tcn)
    if t is None:
        return None
    return CN_CALL_TABLE.get((bool(wgd), t, _norm_cn(lcn)))


__all__ = ["CN_CALL_TABLE", "get_cn_call"]

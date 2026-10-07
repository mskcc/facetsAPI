"""Tiered fit selection for single-run-per-sample mode.

A manifest lists every fit of a sample with a ``review_status``. Fits are grouped in
three tiers, searched in the order the policy gives::

    reviewed : reviewed_best_fit  > reviewed_acceptable_fit        (a human chose it)
    autoqc   : auto_qc_best_fit   > auto_qc_pass                   (the autoQC ruleset chose it)
    default  : fit_name == "default", whatever its status          (facets' own first fit)

``FitPolicy(("reviewed", "autoqc", "default"))`` is the default. Restricting to one
tier is ``FitPolicy.reviewed_only()`` etc. Within a status, manifest order wins unless
``tie_break="newest_reviewed"``.

Selection FALLS THROUGH: the first candidate whose files exist for the requested
purity/hisens mode is taken. Ultra fits (``.../ultra``, hisens-only on disk) are
excluded unless ``include_ultra_fits=True``.
"""

from dataclasses import dataclass

import pandas as pd

from .exceptions import FacetsConfigError
from .layout import NOT_SELECTED, fit_dir, is_ultra, locate_fit_files

TIER_REVIEWED = "reviewed"
TIER_AUTOQC = "autoqc"
TIER_DEFAULT = "default"
TIER_UNREVIEWED = "unreviewed"   # a loaded (all-fits mode) row that belongs to no selection tier

TIERS = {
    TIER_REVIEWED: ("reviewed_best_fit", "reviewed_acceptable_fit"),
    TIER_AUTOQC:   ("auto_qc_best_fit", "auto_qc_pass"),
    TIER_DEFAULT:  None,
}
TIE_BREAKS = ("manifest_order", "newest_reviewed")

STATUS_LABELS = {
    "reviewed_best_fit":       "Best",
    "reviewed_acceptable_fit": "Acceptable",
    "auto_qc_best_fit":        "AutoQC_Best",
    "auto_qc_pass":            "AutoQC_Pass",
}
DEFAULT_LABEL = "Default"
UNREVIEWED_LABEL = "Unreviewed"


def tier_of(review_status, fit_name):
    """Natural tier of a manifest row (independent of any policy)."""
    for tier, statuses in TIERS.items():
        if statuses and review_status in statuses:
            return tier
    if str(fit_name) == "default":
        return TIER_DEFAULT
    return TIER_UNREVIEWED


def label_of(review_status, fit_name):
    """Short label used in fit_map / writeReport: Best, Acceptable, AutoQC_Best, AutoQC_Pass, Default, Unreviewed."""
    if review_status in STATUS_LABELS:
        return STATUS_LABELS[review_status]
    if str(fit_name) == "default":
        return DEFAULT_LABEL
    return UNREVIEWED_LABEL


@dataclass(frozen=True)
class FitCandidate:
    fit_name: str
    review_status: str
    reviewed_by: str
    tier: str
    label: str
    row_index: int

    @property
    def fit_status(self):
        """The legacy 'status#fit_name#reviewed_by' string stored at MetaDictMap.FIT_STATUS."""
        return "%s#%s#%s" % (self.review_status, self.fit_name, self.reviewed_by)


@dataclass(frozen=True)
class FitPolicy:
    tiers: tuple = (TIER_REVIEWED, TIER_AUTOQC, TIER_DEFAULT)
    include_ultra_fits: bool = False
    tie_break: str = "manifest_order"

    def __post_init__(self):
        tiers = tuple(self.tiers) if not isinstance(self.tiers, str) else (self.tiers,)
        object.__setattr__(self, "tiers", tiers)
        if not tiers:
            raise FacetsConfigError("FitPolicy needs at least one tier")
        bad = [t for t in tiers if t not in TIERS]
        if bad:
            raise FacetsConfigError("unknown fit tier(s) %s; choose from %s" % (bad, list(TIERS)))
        if len(set(tiers)) != len(tiers):
            raise FacetsConfigError("duplicate tiers in FitPolicy: %s" % (tiers,))
        if self.tie_break not in TIE_BREAKS:
            raise FacetsConfigError("tie_break must be one of %s" % (TIE_BREAKS,))

    # ---- convenience constructors -------------------------------------------------
    @classmethod
    def reviewed_only(cls, **kw):
        return cls(tiers=(TIER_REVIEWED,), **kw)

    @classmethod
    def autoqc_only(cls, **kw):
        return cls(tiers=(TIER_AUTOQC,), **kw)

    @classmethod
    def default_only(cls, **kw):
        return cls(tiers=(TIER_DEFAULT,), **kw)

    @classmethod
    def reviewed_then_autoqc(cls, **kw):
        return cls(tiers=(TIER_REVIEWED, TIER_AUTOQC), **kw)

    @property
    def allows_defaults(self):
        return TIER_DEFAULT in self.tiers

    def with_defaults(self, allow):
        """Return a copy with the default tier appended or removed (used by setSingleRunPerSample)."""
        tiers = [t for t in self.tiers if t != TIER_DEFAULT]
        if allow:
            tiers.append(TIER_DEFAULT)
        if not tiers:
            tiers = [TIER_DEFAULT]
        return FitPolicy(tuple(tiers), self.include_ultra_fits, self.tie_break)

    # ---- selection ------------------------------------------------------------------
    def eligible_rows(self, manifest_df):
        """Manifest rows that can carry a fit under this policy ('Not selected' and, unless
        include_ultra_fits, ultra rows removed). Keeps the original index."""
        df = manifest_df
        mask = df["fit_name"].astype(str) != NOT_SELECTED
        mask &= df["fit_name"].astype(str) != ""
        if not self.include_ultra_fits:
            mask &= ~df["fit_name"].map(is_ultra)
        return df.loc[mask]

    def _ordered(self, rows):
        if self.tie_break == "newest_reviewed" and "date_reviewed" in rows.columns:
            dates = pd.to_datetime(rows["date_reviewed"], errors="coerce")
            return rows.assign(_d=dates).sort_values("_d", ascending=False, na_position="last", kind="stable")
        return rows

    def candidates(self, manifest_df):
        """Priority-ordered FitCandidate list for a manifest (no file checks)."""
        rows = self.eligible_rows(manifest_df)
        seen = set()
        out = []

        def emit(sub, tier):
            for idx, r in self._ordered(sub).iterrows():
                if idx in seen:
                    continue
                seen.add(idx)
                status = str(r["review_status"])
                fit = str(r["fit_name"])
                out.append(FitCandidate(fit, status, str(r["reviewed_by"]), tier, label_of(status, fit), int(idx)))

        for tier in self.tiers:
            statuses = TIERS[tier]
            if statuses is None:
                emit(rows.loc[rows["fit_name"].astype(str) == "default"], TIER_DEFAULT)
            else:
                for status in statuses:
                    emit(rows.loc[rows["review_status"].astype(str) == status], tier)
        return out

    def select_fit(self, manifest_df, mode, class_dir):
        """Return (candidate, files, tried) for the first candidate whose files exist, else (None, None, tried).

        ``tried`` is the list of candidates whose files were missing (for diagnostics).
        """
        tried = []
        for cand in self.candidates(manifest_df):
            files, _missing = locate_fit_files(fit_dir(class_dir, cand.fit_name), mode)
            if files is not None:
                return cand, files, tried
            tried.append(cand)
        return None, None, tried


__all__ = ["FitPolicy", "FitCandidate", "TIERS", "TIER_REVIEWED", "TIER_AUTOQC", "TIER_DEFAULT",
           "TIER_UNREVIEWED", "STATUS_LABELS", "tier_of", "label_of"]

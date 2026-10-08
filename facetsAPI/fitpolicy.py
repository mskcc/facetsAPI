"""Tiered fit selection for single-run-per-sample mode.

A manifest lists every fit of a sample with a ``review_status``. Fits are grouped in
three tiers, searched in the order the policy gives::

    reviewed : reviewed_best_fit  > reviewed_acceptable_fit  by a HUMAN reviewer
    autoqc   : auto_qc_best_fit   > auto_qc_pass             (current autoQC), and
               reviewed_best_fit  > reviewed_acceptable_fit  whose reviewed_by matches the automated
               reviewer pattern (legacy autoQC wrote 'reviewed_best_fit' as 'auto-qc-script',
               'richara4/auto-qc-script', 'bandlamc/auto-qc-script', ...)
    default  : fit_name == "default", whatever its status          (facets' own first fit)

A row's tier therefore depends on BOTH review_status and reviewed_by.

``FitPolicy(("reviewed", "autoqc", "default"))`` is the default. Restricting to one
tier is ``FitPolicy.reviewed_only()`` etc. Within a status, manifest order wins unless
``tie_break="newest_reviewed"``.

Selection FALLS THROUGH: the first candidate whose files exist for the requested
purity/hisens mode is taken. Ultra fits (``.../ultra``, hisens-only on disk) are
excluded unless ``include_ultra_fits=True``.
"""

import re
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

# reviewed_by values produced by automated review: 'auto-qc', 'auto-qc-script', '<user>/auto-qc-script', 'auto_qc' ...
AUTO_REVIEWER_PATTERN = r"auto[-_]?qc"
REVIEWED_STATUSES = TIERS[TIER_REVIEWED]
AUTOQC_STATUSES = TIERS[TIER_AUTOQC]


def is_auto_reviewer(reviewed_by, pattern=AUTO_REVIEWER_PATTERN):
    return bool(re.search(pattern, str(reviewed_by or ""), flags=re.IGNORECASE))


def classify(review_status, fit_name, reviewed_by=None, pattern=AUTO_REVIEWER_PATTERN):
    """Return (tier, rank, label) for a manifest row. rank orders rows within a tier (0 = best)."""
    status = str(review_status)
    if status in REVIEWED_STATUSES:
        rank = REVIEWED_STATUSES.index(status)
        if is_auto_reviewer(reviewed_by, pattern):
            return TIER_AUTOQC, rank, (STATUS_LABELS["auto_qc_best_fit"], STATUS_LABELS["auto_qc_pass"])[rank]
        return TIER_REVIEWED, rank, STATUS_LABELS[status]
    if status in AUTOQC_STATUSES:
        rank = AUTOQC_STATUSES.index(status)
        return TIER_AUTOQC, rank, STATUS_LABELS[status]
    if str(fit_name) == "default":
        return TIER_DEFAULT, 0, DEFAULT_LABEL
    return TIER_UNREVIEWED, 0, UNREVIEWED_LABEL



def tier_of(review_status, fit_name, reviewed_by=None, pattern=AUTO_REVIEWER_PATTERN):
    """Natural tier of a manifest row (independent of any policy); see classify()."""
    return classify(review_status, fit_name, reviewed_by, pattern)[0]


def label_of(review_status, fit_name, reviewed_by=None, pattern=AUTO_REVIEWER_PATTERN):
    """Short label used in fit_map / writeReport: Best, Acceptable, AutoQC_Best, AutoQC_Pass, Default, Unreviewed."""
    return classify(review_status, fit_name, reviewed_by, pattern)[2]


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
    auto_reviewer_pattern: str = AUTO_REVIEWER_PATTERN   # reviewed_by regex that marks a row as automated

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
        try:
            re.compile(self.auto_reviewer_pattern)
        except re.error as e:
            raise FacetsConfigError("auto_reviewer_pattern is not a valid regex: %s" % e)

    def classify(self, review_status, fit_name, reviewed_by=None):
        return classify(review_status, fit_name, reviewed_by, self.auto_reviewer_pattern)

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
        return FitPolicy(tuple(tiers), self.include_ultra_fits, self.tie_break, self.auto_reviewer_pattern)

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

    def _classified(self, rows):
        """[(index, tier, rank, label)] for eligible rows, in manifest order."""
        out = []
        for idx, r in rows.iterrows():
            tier, rank, label = self.classify(r["review_status"], r["fit_name"], r["reviewed_by"])
            out.append((idx, tier, rank, label))
        return out

    def candidates(self, manifest_df):
        """Priority-ordered FitCandidate list for a manifest (no file checks).

        For each tier in policy order: rows classified into that tier (by status AND reviewer), best rank
        first, ties by tie_break. The default tier takes every fit_name == "default" row not already emitted,
        whatever its status.
        """
        rows = self.eligible_rows(manifest_df)
        info = {idx: (tier, rank, label) for idx, tier, rank, label in self._classified(rows)}
        order = {idx: pos for pos, idx in enumerate(self._ordered(rows).index)}
        seen = set()
        out = []

        def emit(indices):
            for idx in sorted(indices, key=lambda i: (info[i][1], order[i])):
                if idx in seen:
                    continue
                seen.add(idx)
                r = rows.loc[idx]
                tier, _rank, label = info[idx]
                out.append(FitCandidate(str(r["fit_name"]), str(r["review_status"]), str(r["reviewed_by"]),
                                        tier, label, int(idx)))

        for tier in self.tiers:
            if tier == TIER_DEFAULT:
                emit([i for i in rows.index if str(rows.loc[i, "fit_name"]) == "default"])
            else:
                emit([i for i in rows.index if info[i][0] == tier])
        # a default row claimed by the default tier keeps its own label but reports the tier it was found in
        return [FitCandidate(c.fit_name, c.review_status, c.reviewed_by,
                             TIER_DEFAULT if (c.tier != TIER_DEFAULT and c.fit_name == "default" and
                                              info[c.row_index][0] not in self.tiers) else c.tier,
                             c.label, c.row_index) for c in out]

    def representative_rows(self, manifest_df):
        """One row per fit_name from eligible_rows(): manifests keep review history, so a fit can appear
        several times (e.g. 'default' as reviewed_best_fit and again as not_reviewed). The highest-ranked
        row wins: tier (reviewed > autoqc > default > unreviewed, classified by status AND reviewer), then
        rank within the tier, then newest date_reviewed, then the LAST manifest row.
        Returns (rows, n_duplicates_dropped)."""
        rows = self.eligible_rows(manifest_df)
        tier_rank = {TIER_REVIEWED: 0, TIER_AUTOQC: 1, TIER_DEFAULT: 2, TIER_UNREVIEWED: 3}
        dates = pd.to_datetime(rows["date_reviewed"], errors="coerce") if "date_reviewed" in rows.columns else None
        best = {}
        for pos, (idx, tier, rank, _label) in enumerate(self._classified(rows)):
            fit = str(rows.loc[idx, "fit_name"])
            d = dates.loc[idx] if dates is not None else pd.NaT
            key = (tier_rank[tier], rank,
                   0 if pd.isna(d) else -d.value,   # newest review first; undated rows after dated ones
                   -pos)                             # then the later manifest row
            if fit not in best or key < best[fit][0]:
                best[fit] = (key, idx)
        keep = [idx for _k, idx in best.values()]
        return rows.loc[sorted(keep, key=lambda i: rows.index.get_loc(i))], len(rows) - len(keep)

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
           "TIER_UNREVIEWED", "STATUS_LABELS", "AUTO_REVIEWER_PATTERN", "is_auto_reviewer", "classify",
           "tier_of", "label_of"]

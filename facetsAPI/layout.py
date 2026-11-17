"""Repository layout discovery for standard and 2n FACETS samples.

Standard sample directory::

    <repo>/<P-00xxx>/<tumor>_<normal>/facets_review.manifest
                                     /facets_qc.txt
                                     /<fit_name>/...

2n pair directory (no root manifest)::

    <repo>/<P-00xxx>/<tumor>_<normal>/clinical/{facets_review.manifest,facets_qc.txt,<fit>/...}
                                     /research/{facets_review.manifest,facets_qc.txt,<fit>/...}
                                     /research/<fit>/ultra/   (hisens-only fits)

Fit directories are ALWAYS resolved relative to the discovered class directory,
never from the absolute ``path`` column inside the manifest (which records where
facets-preview saw the repository, not where it is mounted now).
"""

import glob
import os
from dataclasses import dataclass

import pandas as pd

from .exceptions import FacetsConfigError, ManifestError

STANDARD = "standard"
CLINICAL = "clinical"
RESEARCH = "research"
TWO_N_CLASSES = (CLINICAL, RESEARCH)
FIT_CLASS_CHOICES = (CLINICAL, RESEARCH, "both")

MANIFEST_NAME = "facets_review.manifest"
FACETS_QC_NAME = "facets_qc.txt"
MANIFEST_REQUIRED = ("path", "review_status", "fit_name", "reviewed_by")
NOT_SELECTED = "Not selected"
MODES = ("purity", "hisens")

# master_file_dict slots that must exist for a fit to be usable.
CRITICAL_FILES = ("out", "cncf", "qc", "gene", "adjseg", "unadjseg")


@dataclass(frozen=True)
class FitClass:
    """One manifest-bearing directory for a sample: the sample root (standard) or a 2n class dir."""
    name: str           # "standard" | "clinical" | "research"
    class_dir: str      # directory holding the manifest (no trailing slash)
    manifest_path: str
    facets_qc_path: str

    @property
    def is_2n(self):
        return self.name != STANDARD


def validate_fit_class(fit_class):
    if fit_class not in FIT_CLASS_CHOICES:
        raise FacetsConfigError("fit_class must be one of %s, got %r" % (", ".join(FIT_CLASS_CHOICES), fit_class))
    return fit_class


def validate_mode(mode):
    if mode not in MODES:
        raise FacetsConfigError("hisens_vs_purity must be 'purity' or 'hisens', got %r" % (mode,))
    return mode


def discover_sample(sample_dir):
    """Return the FitClass list for a sample directory.

    Standard layout -> one FitClass named "standard". 2n layout -> one FitClass per
    class subdirectory that holds a manifest, clinical first. No manifest anywhere -> [].
    """
    sample_dir = os.path.normpath(sample_dir)
    root_manifest = os.path.join(sample_dir, MANIFEST_NAME)
    if os.path.isfile(root_manifest):
        return [FitClass(STANDARD, sample_dir, root_manifest, os.path.join(sample_dir, FACETS_QC_NAME))]
    found = []
    for cls in TWO_N_CLASSES:
        cdir = os.path.join(sample_dir, cls)
        man = os.path.join(cdir, MANIFEST_NAME)
        if os.path.isfile(man):
            found.append(FitClass(cls, cdir, man, os.path.join(cdir, FACETS_QC_NAME)))
    return found


def select_classes(classes, fit_class):
    """Apply the user's fit_class choice. Standard samples are unaffected."""
    if not classes or not classes[0].is_2n or fit_class == "both":
        return list(classes)
    return [c for c in classes if c.name == fit_class]


def count_leading_comment_lines(path):
    n = 0
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                n += 1
            else:
                break
    return n


def read_manifest(path):
    """Read a facets_review.manifest, tolerating an optional leading '#' comment block.

    Returns a DataFrame with at least MANIFEST_REQUIRED columns; ``reviewed_by`` NAs are
    normalized to the string "NA" so fit-status strings are stable.
    """
    try:
        skip = count_leading_comment_lines(path)
        df = pd.read_csv(path, sep="\t", skiprows=skip, low_memory=False, dtype=str, keep_default_na=False,
                         na_values=[""])
    except OSError as e:
        raise ManifestError("cannot read manifest %s: %s" % (path, e)) from e
    except (pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        raise ManifestError("malformed manifest %s: %s" % (path, e)) from e
    missing = [c for c in MANIFEST_REQUIRED if c not in df.columns]
    if missing:
        raise ManifestError("manifest %s lacks required column(s) %s (found: %s)"
                            % (path, ", ".join(missing), ", ".join(map(str, df.columns))))
    df = df.copy()
    df["fit_name"] = df["fit_name"].fillna("").astype(str).str.strip()
    df["review_status"] = df["review_status"].fillna("").astype(str).str.strip()
    df["reviewed_by"] = df["reviewed_by"].where(df["reviewed_by"].notna(), "NA").astype(str)
    return df


def short_sample_id(key):
    """Reduce any master_file_dict key or run id to the tumor sample id.

    '<tumor>' | '<tumor>#<class>' | '<tumor>_<normal>' | '<tumor>_<normal>#<class>/<fit>' -> '<tumor>'
    """
    return str(key).split("#")[0].split("_")[0]


def pair_id(key):
    """'<tumor>_<normal>#...' -> '<tumor>_<normal>'; a bare tumor id is returned unchanged."""
    return str(key).split("#")[0]


def is_ultra(fit_name):
    name = str(fit_name).rstrip("/")
    return name == "ultra" or name.endswith("/ultra")


def fit_dir(class_dir, fit_name):
    """The on-disk directory for a fit, resolved from the DISCOVERED class directory."""
    return os.path.join(class_dir, str(fit_name).strip("/"))


def _first(paths):
    paths = sorted(paths)
    return paths[0] if paths else ""


def locate_fit_files(fit_directory, mode):
    """Find the files a fit needs for ``mode`` ('purity' or 'hisens').

    Returns (files, missing). ``files`` is a dict with keys out, cncf, adjseg, unadjseg, rdata, qc,
    gene, ccf, nonsignedout (optional ones are "" when absent); ``missing`` lists the critical
    slots that were not found. ``files`` is None when anything critical is missing.
    """
    validate_mode(mode)
    d = fit_directory.rstrip("/")
    files = {
        "out":      _first(glob.glob(os.path.join(d, "*_%s.out" % mode))),
        "cncf":     _first(glob.glob(os.path.join(d, "*_%s.cncf.txt" % mode))),
        "adjseg":   _first(glob.glob(os.path.join(d, "*_%s_diplogR.adjusted.seg" % mode))),
        "unadjseg": _first(glob.glob(os.path.join(d, "*_%s_diplogR.unadjusted.seg" % mode))),
        "rdata":    _first(glob.glob(os.path.join(d, "*_%s.Rdata" % mode))),
        "qc":       _first(glob.glob(os.path.join(d, "*.qc.txt"))),
        "gene":     _first(glob.glob(os.path.join(d, "*gene_level.txt"))),
        "ccf":      "",
        "nonsignedout": "",
    }
    for maf in sorted(glob.glob(os.path.join(d, "*ccf.maf"))):
        if "nonsignedout" in os.path.basename(maf):
            files["nonsignedout"] = maf
        else:
            files["ccf"] = maf
    missing = [k for k in CRITICAL_FILES if not files[k]]
    return (None if missing else files), missing


__all__ = [
    "FitClass", "discover_sample", "select_classes", "read_manifest", "short_sample_id", "pair_id",
    "is_ultra", "fit_dir", "locate_fit_files", "validate_fit_class", "validate_mode",
    "STANDARD", "CLINICAL", "RESEARCH", "TWO_N_CLASSES", "FIT_CLASS_CHOICES", "NOT_SELECTED",
    "MANIFEST_REQUIRED", "CRITICAL_FILES", "MODES",
]

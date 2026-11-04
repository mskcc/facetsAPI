"""JSON request spec for facetsAPI jobs, and the script generator behind the FacetFlow dashboard.

A spec describes one job: where the samples come from, which fits and filters apply, and which
outputs to write. validate() checks it and render() turns it into a runnable Python script.
See doc/spec.md for the format.

Only the standard library is used and there are no package-relative imports, so a service can load
this file by path without pandas. Rendered scripts run code from this file and nothing else: every
value from a spec is emitted as a checked literal, so a spec cannot inject code.

Hard rule: a job reads FACETS repositories and never writes, modifies or deletes anything in them. No output
type here may call a method that writes to a repository, and every script starts with facetsAPI.write_guard,
which enforces this while the job runs.
"""
import argparse
import json
import math
import re
import sys

SPEC_VERSION = 1

MODES = ("purity", "hisens")
FIT_CLASSES = ("clinical", "research", "both")
TIERS = ("reviewed", "autoqc", "default")
TIE_BREAKS = ("manifest_order", "newest_reviewed")
TRISTATE_FILTERS = {"wgd": "setWGDFilter", "facets_qc": "setFacetsQCFilter"}

#key: (setter, lowest allowed, highest allowed). A missing bound means no limit on that side.
RANGE_FILTERS = {
    "purity":          ("setPurityFilter", 0.0, 1.0),
    "clinical_purity": ("setClinicalPurityFilter", 0.0, 100.0),
    "ploidy":          ("setPloidyFilter", None, None),
    "diplogr":         ("setDipLogRFilter", None, None),
    "cval":            ("setCvalFilter", None, None),
    "fga":             ("setFGAFilter", 0.0, 1.0),
    "frac_loh":        ("setFracLohFilter", 0.0, 1.0),
    "tmb":             ("setTMBFilter", None, None),
    "msi":             ("setMSIFilter", None, None),
}
LIST_FILTERS = {"cancer_type": "setCancerTypeFilter", "cancer_type_detailed": "setCancerTypeDetailFilter"}

#Attribute names on FacetsMeta holding lists of OncoTree codes.
ONCOTREE_PRESETS = {
    "breast_carcinoma": "Breast carcinoma",
    "nscLung_cancer": "Non-small cell lung cancer",
    "endometrial_cancer": "Endometrial cancer",
    "colorectal_cancer": "Colorectal cancer",
    "melanoma": "Melanoma",
    "scLung_cancer": "Small cell lung cancer",
    "ovarian_cancer": "Ovarian cancer",
    "prostate_cancer": "Prostate cancer",
    "unknown_primary": "Cancer of unknown primary",
    "bladder_cancer": "Bladder cancer",
    "pancreatic_cancer": "Pancreatic cancer",
    "renalCell_carcinoma": "Renal cell carcinoma",
    "glioma": "Glioma",
    "headNeck_carcinoma": "Head and neck carcinoma",
    "germCell_tumor": "Germ cell tumor",
}

HISTOGRAM_VARIABLES = ("purity", "clinicalPurity", "ploidy", "dipLogR", "cval", "fga", "frac_loh", "tmb", "msi",
                       "cancerType", "onkoCode", "wgd", "review_status", "fit_class", "fit_name")

#key: (MetaDictMap attribute, output file name)
MERGED_FILES = {
    "adjusted_seg":   ("ADJUSTED_SEG_FILE", "merged_adjusted.seg"),
    "unadjusted_seg": ("UNADJUSTED_SEG_FILE", "merged_unadjusted.seg"),
    "cncf":           ("CNCF_FILE", "merged_cncf.txt"),
    "gene_level":     ("GENE_FILE", "merged_gene_level.txt"),
    "qc":             ("QC_FILE", "merged_qc.txt"),
    "facets_qc":      ("FACETS_QC_FILE", "merged_facets_qc.txt"),
    "ccf_maf":        ("CCF_MAF_FILE", "merged_ccf.maf"),
    "nonsigned_maf":  ("NONSIGNED_MAF_FILE", "merged_nonsigned.maf"),
}

#key: (class holding the setting, attribute, value in facetsAPI)
ALTERATION_THRESHOLDS = {
    "alteration_arm_num": ("FacetsRun", "alteration_arm_num", 20),
    "alteration_cov_min": ("FacetsRun", "alteration_cov_min", 50),
    "percent_arm_loh":    ("FacetsSegment", "percent_arm_loh", 50.0),
    "min_gain_tcn":       ("FacetsSegment", "min_gain_tcn", 4),
    "percent_arm_gain":   ("FacetsSegment", "percent_arm_gain", 50.0),
}

OUTPUT_LABELS = {"cin": "CIN scores", "chromothripsis": "Chromothripsis", "state_percents": "HETLOSS and HOMDEL percentages",
                 "ascets": "ASCETS"}
OUTPUT_TYPES = ("report", "gene_level", "alteration", "summary", "histograms", "purity_cf", "merged",
                "cin", "chromothripsis", "state_percents", "ascets", "copy_folders")

#Each sample must appear once in these: their results are keyed by sample.
ONE_FIT_OUTPUTS = ("cin", "chromothripsis", "state_percents", "ascets")
#These read run.segments, which runAlterationAnalysis() splits at the centromeres, so they are written first.
SEGMENT_SCORE_OUTPUTS = ("cin", "chromothripsis", "state_percents")

CIN_LEVELS = ("full", "chrom", "arm")
GENE_OR_SEG = ("seg", "gene")
BG_METRICS = ("median", "mean")
ASCETS_GENOMES = ("hg19", "hg38")

MAX_SAMPLE_IDS = 100000
LARGE_SAMPLE_LIST = 500

_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
_SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,99}$")
_PATIENT_RE = re.compile(r"^P-\d{7}$")
_GENE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,49}$")


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _is_choice(value, choices):
    return isinstance(value, str) and value in choices


def _is_clean_text(value):
    return isinstance(value, str) and bool(value.strip()) and value.isprintable()


def _is_clean_path(value):
    return _is_clean_text(value) and value.startswith("/")


def _section(spec, key):
    value = spec.get(key)
    return value if isinstance(value, dict) else {}


def _outputs(spec):
    outputs = spec.get("outputs")
    return [o for o in outputs if isinstance(o, dict)] if isinstance(outputs, list) else []


def _has_clinical_filter(filters):
    keys = ("clinical_purity", "tmb", "msi", "cancer_type", "cancer_type_detailed", "oncotree")
    return any(key in filters for key in keys)


def validate(spec):
    """Check a spec. Returns (errors, warnings), both lists of strings. A spec with errors must not be rendered."""
    if not isinstance(spec, dict):
        return ["The spec must be an object."], []
    errors, warnings = [], []

    if spec.get("spec_version") != SPEC_VERSION:
        errors.append("spec_version must be %d." % SPEC_VERSION)
    name = spec.get("name")
    if not (isinstance(name, str) and _NAME_RE.match(name)):
        errors.append("Job name must be 1-64 letters, numbers, dots, dashes or underscores, starting with a letter or number.")

    _check_source(spec, errors)
    _check_fits(spec, errors)
    _check_samples(spec, errors)
    _check_filters(spec, errors, warnings)
    _check_outputs(spec, errors)

    if not errors:
        warnings.extend(_scope_warnings(spec))
    return errors, warnings


def _check_source(spec, errors):
    source = spec.get("source")
    if not isinstance(source, dict):
        errors.append("source is required.")
        return
    if not _is_clean_path(source.get("facets_dir")):
        errors.append("source.facets_dir must be an absolute path.")
    clinical = source.get("clinical_file", "")
    if clinical != "" and not _is_clean_path(clinical):
        errors.append("source.clinical_file must be an absolute path or empty.")
    if clinical == "" and _has_clinical_filter(_section(spec, "filters")):
        errors.append("Cancer type, oncotree, clinical purity, TMB and MSI filters need a clinical file.")
    if source.get("mode", "purity") not in MODES:
        errors.append("source.mode must be one of: %s." % ", ".join(MODES))
    if source.get("fit_class", "research") not in FIT_CLASSES:
        errors.append("source.fit_class must be one of: %s." % ", ".join(FIT_CLASSES))


def _check_fits(spec, errors):
    fits = spec.get("fits")
    if not isinstance(fits, dict):
        errors.append("fits is required.")
        return
    if not isinstance(fits.get("single_run"), bool):
        errors.append("fits.single_run must be true or false.")
    tiers = fits.get("tiers")
    if not (isinstance(tiers, list) and tiers and all(isinstance(t, str) and t in TIERS for t in tiers)
            and len(set(tiers)) == len(tiers)):
        errors.append("fits.tiers must be a non-empty list without repeats, from: %s." % ", ".join(TIERS))
    if not _is_choice(fits.get("tie_break", TIE_BREAKS[0]), TIE_BREAKS):
        errors.append("fits.tie_break must be one of: %s." % ", ".join(TIE_BREAKS))


def _check_samples(spec, errors):
    samples = spec.get("samples")
    if not isinstance(samples, dict):
        errors.append("samples is required.")
        return
    source = samples.get("source")
    if source == "all":
        return
    if source != "list":
        errors.append("samples.source must be 'all' or 'list'.")
        return
    ids = samples.get("ids")
    if not (isinstance(ids, list) and ids):
        errors.append("Enter at least one sample ID, or select all samples.")
        return
    if len(ids) > MAX_SAMPLE_IDS:
        errors.append("At most %d sample IDs per job." % MAX_SAMPLE_IDS)
        return
    patients = [i for i in ids if isinstance(i, str) and _PATIENT_RE.match(i)]
    if patients:
        errors.append("Patient IDs don't match any sample (%d found, e.g. %s). Use tumor sample IDs or pair IDs."
                      % (len(patients), patients[0]))
    invalid = [i for i in ids if not (isinstance(i, str) and _SAMPLE_RE.match(i))]
    if invalid:
        errors.append("%d sample IDs contain invalid characters." % len(invalid))


def _check_filters(spec, errors, warnings):
    filters = spec.get("filters", {})
    if not isinstance(filters, dict):
        errors.append("filters must be an object.")
        return
    known = set(RANGE_FILTERS) | set(TRISTATE_FILTERS) | set(LIST_FILTERS) | {"oncotree"}
    for key in filters:
        if key not in known:
            errors.append("Unknown filter: %s." % key)
    for key in RANGE_FILTERS:
        if key in filters:
            _check_range(key, filters[key], errors)
    for key in TRISTATE_FILTERS:
        if key in filters and not isinstance(filters[key], bool):
            errors.append("Filter %s must be true or false." % key)
    for key in LIST_FILTERS:
        if key in filters and not _is_text_list(filters[key]):
            errors.append("Filter %s must be a non-empty list of values." % key)
    if "oncotree" in filters:
        _check_oncotree(filters["oncotree"], errors)

    clinical = filters.get("clinical_purity")
    if isinstance(clinical, dict) and _is_number(clinical.get("max")) and clinical["max"] <= 1:
        warnings.append("Clinical purity is a percentage (0-100). A maximum of %s keeps almost no samples." % clinical["max"])


def _is_text_list(value):
    return isinstance(value, list) and bool(value) and all(_is_clean_text(v) for v in value)


def _check_range(key, bounds, errors):
    if not isinstance(bounds, dict) or not ("min" in bounds or "max" in bounds):
        errors.append("Filter %s needs a minimum or a maximum." % key)
        return
    lowest, highest = RANGE_FILTERS[key][1:]
    for side in ("min", "max"):
        if side not in bounds:
            continue
        value = bounds[side]
        if not _is_number(value):
            errors.append("Filter %s: %s must be a number." % (key, side))
        elif (lowest is not None and value < lowest) or (highest is not None and value > highest):
            errors.append("Filter %s: %s must be between %s and %s." % (key, side, lowest, highest))
    low, high = bounds.get("min"), bounds.get("max")
    if _is_number(low) and _is_number(high) and low > high:
        errors.append("Filter %s: minimum is above the maximum." % key)


def _check_oncotree(oncotree, errors):
    if not isinstance(oncotree, dict):
        errors.append("Filter oncotree must be an object.")
        return
    presets = oncotree.get("presets", [])
    codes = oncotree.get("codes", [])
    if not (isinstance(presets, list) and all(isinstance(p, str) and p in ONCOTREE_PRESETS for p in presets)):
        errors.append("Filter oncotree: unknown preset.")
    if not (isinstance(codes, list) and all(_is_clean_text(c) for c in codes)):
        errors.append("Filter oncotree: codes must be a list of values.")
    if not presets and not codes:
        errors.append("Filter oncotree needs at least one preset or code.")


def _check_outputs(spec, errors):
    outputs = spec.get("outputs")
    if not (isinstance(outputs, list) and outputs):
        errors.append("Choose at least one output.")
        return
    seen = set()
    for output in outputs:
        kind = output.get("type") if isinstance(output, dict) else None
        if kind not in OUTPUT_TYPES:
            errors.append("Unknown output type: %s." % (kind,))
            continue
        if kind in seen:
            errors.append("Output %s is listed twice." % kind)
        seen.add(kind)
        if kind == "gene_level":
            genes = output.get("genes", [])
            if not (isinstance(genes, list) and all(isinstance(g, str) and _GENE_RE.match(g) for g in genes)):
                errors.append("Gene-level report: gene names may only contain letters, numbers, dots, dashes and underscores.")
        elif kind == "alteration":
            _check_alteration(output, errors)
        elif kind == "histograms":
            variables = output.get("variables")
            if not (isinstance(variables, list) and variables and all(_is_choice(v, HISTOGRAM_VARIABLES) for v in variables)):
                errors.append("Histograms: choose at least one variable from: %s." % ", ".join(HISTOGRAM_VARIABLES))
            bins = output.get("bins", 10)
            if not (isinstance(bins, int) and not isinstance(bins, bool) and 1 <= bins <= 200):
                errors.append("Histograms: bins must be a whole number from 1 to 200.")
        elif kind == "merged":
            files = output.get("files")
            if not (isinstance(files, list) and files and all(_is_choice(f, MERGED_FILES) for f in files)):
                errors.append("Merged files: choose at least one from: %s." % ", ".join(MERGED_FILES))
        elif kind == "cin":
            levels = output.get("levels")
            if not (isinstance(levels, list) and levels and all(_is_choice(v, CIN_LEVELS) for v in levels)
                    and len(set(levels)) == len(levels)):
                errors.append("CIN scores: choose at least one level from: %s." % ", ".join(CIN_LEVELS))
            if not _is_choice(output.get("gene_or_seg", "seg"), GENE_OR_SEG):
                errors.append("CIN scores: calculate from segments or genes.")
            if not isinstance(output.get("split_arms", True), bool):
                errors.append("CIN scores: split_arms must be true or false.")
        elif kind == "chromothripsis":
            _check_chromothripsis(output, errors)
        elif kind == "state_percents":
            if not (isinstance(output.get("include_x", False), bool) and isinstance(output.get("include_y", False), bool)):
                errors.append("HETLOSS and HOMDEL percentages: include_x and include_y must be true or false.")
        elif kind == "ascets":
            _check_ascets(output, errors)

    one_fit = [k for k in ONE_FIT_OUTPUTS if k in seen]
    if one_fit:
        fits = _section(spec, "fits")
        if fits.get("single_run") is False or _section(spec, "source").get("fit_class") == "both":
            errors.append("%s need one fit per sample: turn on \"One fit per sample\" and don't use both 2n fit classes."
                          % ", ".join(OUTPUT_LABELS[k] for k in one_fit))


def _check_chromothripsis(output, errors):
    if not _is_choice(output.get("gene_or_seg", "seg"), GENE_OR_SEG):
        errors.append("Chromothripsis: calculate from segments or genes.")
    if not _is_choice(output.get("bg_metric", "median"), BG_METRICS):
        errors.append("Chromothripsis: background must be median or mean.")
    for key, whole in (("min_breakpoints", True), ("breakpoint_offset", True), ("min_ris_delta_pct", False)):
        value = output.get(key, 0)
        if not _is_number(value) or value < 0 or (whole and value != int(value)):
            errors.append("Chromothripsis: %s must be %s." % (key, "a whole number, 0 or more" if whole else "a number, 0 or more"))


def _check_ascets(output, errors):
    if not _is_choice(output.get("genome", "hg19"), ASCETS_GENOMES):
        errors.append("ASCETS: genome must be hg19 or hg38.")
    for key in ("min_arm_breadth", "arm_alt_frac_thresh"):
        value = output.get(key, 0.5)
        if not _is_number(value) or not 0 <= value <= 1:
            errors.append("ASCETS: %s must be between 0 and 1." % key)
    for key in ("keep_noise", "use_adjusted"):
        if not isinstance(output.get(key, False), bool):
            errors.append("ASCETS: %s must be true or false." % key)


def _check_alteration(output, errors):
    thresholds = output.get("thresholds", {})
    if not isinstance(thresholds, dict):
        errors.append("Alteration analysis: thresholds must be an object.")
        return
    for key, value in thresholds.items():
        if key not in ALTERATION_THRESHOLDS:
            errors.append("Alteration analysis: unknown threshold %s." % key)
        elif not _is_number(value):
            errors.append("Alteration analysis: %s must be a number." % key)


def _scope_warnings(spec):
    warnings = []
    whole_repo = spec["samples"]["source"] == "all"
    outputs = {o["type"]: o for o in _outputs(spec)}
    if whole_repo and not spec["fits"]["single_run"]:
        warnings.append("Whole repository with all fits loads every fit of every sample. This can take many hours "
                        "and a lot of memory; consider one fit per sample.")
    if whole_repo and "gene_level" in outputs and not outputs["gene_level"].get("genes"):
        warnings.append("A gene-level report for all genes over the whole repository produces a very large file; "
                        "consider listing genes.")
    if whole_repo and "copy_folders" in outputs:
        warnings.append("Copying sample folders for the whole repository needs a very large amount of disk space.")
    return warnings


def suggest_resources(spec):
    """Suggest a resource preset for a valid spec: 'small', 'medium' or 'large'."""
    samples = spec["samples"]
    if samples["source"] == "list":
        return "small" if len(samples["ids"]) <= LARGE_SAMPLE_LIST else "medium"
    outputs = _outputs(spec)
    genes_everywhere = any(o["type"] == "gene_level" and not o.get("genes") for o in outputs)
    copies = any(o["type"] == "copy_folders" for o in outputs)
    return "large" if genes_everywhere or copies or not spec["fits"]["single_run"] else "medium"


def samples_text(spec):
    """Contents of samples.txt for a valid spec with a sample list, else None."""
    samples = spec["samples"]
    return "\n".join(samples["ids"]) + "\n" if samples["source"] == "list" else None


def _comment(value):
    return "".join(ch if ch.isprintable() else "?" for ch in str(value))


def _fit_policy_code(fits):
    tiers = fits["tiers"]
    tie_break = fits.get("tie_break", TIE_BREAKS[0])
    extra = "" if tie_break == TIE_BREAKS[0] else "tie_break=%r" % tie_break
    named = {("reviewed",): "reviewed_only", ("autoqc",): "autoqc_only", ("default",): "default_only"}.get(tuple(tiers))
    if named:
        return "FitPolicy.%s(%s)" % (named, extra)
    return "FitPolicy(%r%s)" % (tuple(tiers), ", " + extra if extra else "")


def _bound(value, fallback):
    return repr(value) if value is not None else fallback


def _filter_lines(filters):
    lines = []
    for key, (setter, lowest, highest) in RANGE_FILTERS.items():
        if key in filters:
            low = _bound(filters[key].get("min"), repr(lowest) if lowest is not None else 'float("-inf")')
            high = _bound(filters[key].get("max"), repr(highest) if highest is not None else 'float("inf")')
            lines.append("dataset.%s(%s, %s)" % (setter, low, high))
    for key, setter in TRISTATE_FILTERS.items():
        if key in filters:
            lines.append("dataset.%s(%r)" % (setter, filters[key]))
    for key, setter in LIST_FILTERS.items():
        if key in filters:
            lines.append("dataset.%s(%r)" % (setter, filters[key]))
    if "oncotree" in filters:
        parts = ["FacetsMeta.%s" % p for p in filters["oncotree"].get("presets", [])]
        if filters["oncotree"].get("codes"):
            parts.append(repr(filters["oncotree"]["codes"]))
        lines.append("dataset.setOnkoCodeFilter(%s)" % " + ".join(parts))
    return lines


def _output_lines(output, ctx):
    kind = output["type"]
    if kind == "cin":
        return ["write_cin_results(cin_for_runs(dataset.runList, method=\"Eckert\", gene_or_seg=%r, splitArms=%r, "
                "cinLevel=%r), os.path.join(OUT, %r))"
                % (output.get("gene_or_seg", "seg"), output.get("split_arms", True), level, "cin_%s.tsv" % level)
                for level in output["levels"]]
    if kind == "chromothripsis":
        return ['chromothripsis_for_runs(dataset.runList, os.path.join(OUT, "chromothripsis"), gene_or_seg=%r, '
                "min_breakpoints=%r, min_ris_delta_pct=%r, breakpoint_offset=%r, bg_metric=%r)"
                % (output.get("gene_or_seg", "seg"), output.get("min_breakpoints", 2), output.get("min_ris_delta_pct", 0.25),
                   output.get("breakpoint_offset", 2), output.get("bg_metric", "median"))]
    if kind == "state_percents":
        return ['write_summary(os.path.join(OUT, "hetloss_homdel_percent.tsv"), [{"sample": run.id, '
                "**compute_sample_percentages(run, include_x=%r, include_y=%r)} for run in dataset.runList])"
                % (output.get("include_x", False), output.get("include_y", False))]
    if kind == "ascets":
        return ['ExtTools(meta, config=FacetsConfig.from_env(r_prefix=%r), dataset=dataset).runAscets(os.path.join(OUT, "ascets", ""), '
                "ref_genome_coords=%r, min_arm_breadth=%r, keep_noise=%r, arm_alt_frac_thresh=%r, use_adjusted=%r)"
                % (ctx.get("r_prefix", ""), output.get("genome", "hg19"), output.get("min_arm_breadth", 0.5),
                   output.get("keep_noise", False), output.get("arm_alt_frac_thresh", 0.7), output.get("use_adjusted", True))]
    if kind == "copy_folders":
        return ['dataset.copyDatasetToFolder(os.path.join(OUT, "facets_folders"))']
    if kind == "report":
        return ['dataset.writeReport(os.path.join(OUT, "report.txt"))']
    if kind == "gene_level":
        genes = output.get("genes") or None
        return ['dataset.writeGeneLevelReport(os.path.join(OUT, "gene_level_report.txt"), %r)' % (genes,)]
    if kind == "alteration":
        lines = []
        for key, value in output.get("thresholds", {}).items():
            owner, attribute = ALTERATION_THRESHOLDS[key][:2]
            lines.append("%s.%s = %r" % (owner, attribute, value))
        lines.append("dataset.runAlterationAnalysis()")
        lines.append('dataset.writeAlterationData(os.path.join(OUT, "alteration_data.txt"))')
        return lines
    if kind == "summary":
        return ["dataset.writeDatasetSummary(OUT)"]
    if kind == "histograms":
        return ["dataset.createHistogram(%r, %r, OUT)" % (variable, output.get("bins", 10))
                for variable in output["variables"]]
    if kind == "purity_cf":
        return ['dataset.writePurityCFs(os.path.join(OUT, "purity_cf.txt"))']
    #merged
    lines = ["tools = ExtTools(meta, dataset=dataset)"]
    for key in output["files"]:
        attribute, filename = MERGED_FILES[key]
        lines.append("tools.makeMergedFile(MetaDictMap.%s, os.path.join(OUT, %r))" % (attribute, filename))
    return lines


def render(spec, ctx):
    """Turn a valid spec into a Python script.

    ctx holds: job_id, user, facetsapi_home, facetsapi_sha, output_dir, job_dir, protected_roots (folders no job
    may ever change), and r_prefix when ASCETS is used.
    Raises ValueError if the spec has validation errors.
    """
    errors, _ = validate(spec)
    if errors:
        raise ValueError("Invalid spec: " + " ".join(errors))

    source = spec["source"]
    fits = spec["fits"]
    filters = spec.get("filters", {})
    outputs = spec["outputs"]
    kinds = {o["type"] for o in outputs}

    imports = ["FacetsDataset", "FacetsMeta"]
    if fits["single_run"]:
        imports.append("FitPolicy")
    if "merged" in kinds:
        imports += ["ExtTools", "MetaDictMap"]
    if "ascets" in kinds:
        imports += ["ExtTools", "FacetsConfig"]
    for output in outputs:
        if output["type"] == "alteration":
            imports += sorted({ALTERATION_THRESHOLDS[k][0] for k in output.get("thresholds", {})})

    lines = [
        "# Generated by the FacetFlow FacetsAPI dashboard. Do not edit; it is regenerated from spec.json.",
        "# job: %s" % _comment(ctx["job_id"]),
        "# name: %s" % _comment(spec["name"]),
        "# user: %s" % _comment(ctx["user"]),
        "# spec_version: %d" % SPEC_VERSION,
        "# facetsAPI: %s" % _comment(ctx["facetsapi_sha"]),
        "import json",
        "import os",
        "import sys",
        "",
        "sys.dont_write_bytecode = True",
        "sys.path.insert(0, %r)" % ctx["facetsapi_home"],
        "",
        "OUT = %r" % ctx["output_dir"],
        "JOB_DIR = %r" % ctx["job_dir"],
        'os.environ["MPLCONFIGDIR"] = os.path.join(JOB_DIR, ".matplotlib")',
        "",
        "# Hard rule: this job never writes, modifies or deletes anything in a FACETS repository.",
        "from facetsAPI.write_guard import enable_write_guard",
        "enable_write_guard([OUT, JOB_DIR], protected=%r, allow_ascets=%r)" % (list(ctx.get("protected_roots", [])), "ascets" in kinds),
        "os.makedirs(OUT, exist_ok=True)",
        "",
        "from facetsAPI import %s" % ", ".join(sorted(set(imports))),
    ]
    script_imports = []
    if "cin" in kinds:
        script_imports.append("from cin_scores import cin_for_runs, write_cin_results")
    if "chromothripsis" in kinds:
        script_imports.append("from chromothripsis_scores import chromothripsis_for_runs")
    if "state_percents" in kinds:
        script_imports.append("from get_state_percents import compute_sample_percentages, write_summary")
    if script_imports:
        lines.append("sys.path.insert(0, %r)" % (ctx["facetsapi_home"].rstrip("/") + "/scripts"))
        lines += script_imports
    lines += [
        "",
        "meta = FacetsMeta(%r, %r, %r, fit_class=%r)" % (
            source.get("clinical_file", ""), source["facets_dir"], source.get("mode", "purity"),
            source.get("fit_class", "research")),
    ]
    if fits["single_run"]:
        lines.append("meta.setSingleRunPerSample(True, %r)" % ("default" in fits["tiers"]))
        lines.append("meta.setFitPolicy(%s)" % _fit_policy_code(fits))
    if spec["samples"]["source"] == "list":
        lines.append('meta.selectSamplesFromFile(os.path.join(JOB_DIR, "samples.txt"))')
    lines += ["meta.buildFacetsMeta()", "", "dataset = FacetsDataset(meta)"]
    lines += _filter_lines(filters)
    lines += ["dataset.buildFacetsDataset()", ""]
    for output in sorted(outputs, key=lambda o: o["type"] not in SEGMENT_SCORE_OUTPUTS):
        lines += _output_lines(output, ctx)
    lines += [
        "",
        "summary = {",
        '    "num_runs": len(dataset.runList),',
        '    "num_filtered": dataset.num_filtered,',
        '    "num_file_failed": dataset.num_file_failed,',
        '    "num_build_failed": dataset.num_build_failed,',
        '    "build_counts": meta.build_counts,',
        '    "unmatched_sample_ids": list(meta.unmatched_sample_ids),',
        '    "outputs": sorted(os.listdir(OUT)),',
        "}",
        'with open(os.path.join(JOB_DIR, "run_summary.json"), "w") as fh:',
        "    json.dump(summary, fh, indent=2, default=str)",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate a facetsAPI request spec or render it to a script.")
    parser.add_argument("command", choices=("validate", "render"))
    parser.add_argument("spec", help="Path to the spec JSON file.")
    parser.add_argument("--output-dir", default="./facets_output", help="Where the rendered script writes outputs.")
    parser.add_argument("--job-dir", default=".", help="Folder holding samples.txt; the script writes run_summary.json here.")
    parser.add_argument("--facetsapi-home", default=".", help="facetsAPI checkout the rendered script imports.")
    args = parser.parse_args(argv)

    with open(args.spec) as fh:
        spec = json.load(fh)
    errors, warnings = validate(spec)
    for warning in warnings:
        print("warning: " + warning, file=sys.stderr)
    for error in errors:
        print("error: " + error, file=sys.stderr)
    if errors:
        return 1
    if args.command == "render":
        print(render(spec, {"job_id": "cli", "user": "cli", "facetsapi_home": args.facetsapi_home,
                            "facetsapi_sha": "unknown", "output_dir": args.output_dir, "job_dir": args.job_dir,
                            "protected_roots": [spec["source"]["facets_dir"]]}), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())

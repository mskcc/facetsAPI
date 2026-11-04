# Request spec

A request spec is a JSON description of one facetsAPI job. `facetsAPI/spec.py` validates it and renders it to a Python script. The FacetFlow FacetsAPI dashboard builds these specs from its form, and stores them so a job can be rerun or edited later.

```json
{
  "spec_version": 1,
  "name": "tp53_high_purity",
  "source": {
    "facets_dir": "/path/to/impact_2n/facets/all",
    "clinical_file": "/path/to/data_clinical_sample.oncokb.txt",
    "mode": "purity",
    "fit_class": "research"
  },
  "fits": {"single_run": true, "tiers": ["reviewed"]},
  "samples": {"source": "all"},
  "filters": {"purity": {"min": 0.3}},
  "outputs": [{"type": "report"}, {"type": "gene_level", "genes": ["TP53"]}]
}
```

### Fields

| Field | Values |
|---|---|
| `name` | 1-64 letters, numbers, dots, dashes or underscores |
| `source.mode` | `purity` (default), `hisens` |
| `source.fit_class` | `clinical`, `research` (default), `both`; only affects 2n pairs |
| `source.clinical_file` | Path, or `""` to walk the repository. Clinical filters need a file. |
| `fits.single_run` | `true`: one fit per sample, chosen by `tiers`. `false`: every fit. |
| `fits.tiers` | Ordered subset of `reviewed`, `autoqc`, `default`. `["reviewed"]` means human-reviewed fits only. |
| `fits.tie_break` | `manifest_order` (default), `newest_reviewed`: which fit wins when several in a tier qualify |
| `samples` | `{"source": "all"}`, or `{"source": "list", "ids": [...]}` with tumor sample IDs or pair IDs. Patient IDs are rejected. |

### Filters

All filters are optional. Each one you include must pass (AND).

| Filter | Form | Notes |
|---|---|---|
| `purity`, `fga`, `frac_loh` | `{"min": x, "max": y}` | Fractions, 0-1. `purity` always excludes NA and exactly 0.3, the value FACETS reports when it cannot estimate purity. |
| `clinical_purity` | `{"min": x, "max": y}` | **Percent, 0-100.** |
| `ploidy`, `diplogr`, `cval`, `tmb`, `msi` | `{"min": x, "max": y}` | |
| `wgd`, `facets_qc` | `true` or `false` | Leave out to keep both. |
| `cancer_type`, `cancer_type_detailed` | List of values | Exact, case-sensitive match. |
| `oncotree` | `{"presets": [...], "codes": [...]}` | Presets are the code lists on `FacetsMeta`, such as `breast_carcinoma`. |

Range filters need at least one bound and include the bound itself. A missing bound means no limit on that side.

### Outputs

At least one is required, and each type may appear once.

| Type | Options | Writes |
|---|---|---|
| `report` | | `report.txt` |
| `gene_level` | `genes`: list; empty means all genes | `gene_level_report.txt` |
| `alteration` | `thresholds`: `alteration_arm_num`, `alteration_cov_min`, `percent_arm_loh`, `min_gain_tcn`, `percent_arm_gain` | `alteration_data.txt` |
| `summary` | | `FacetsDataset_summary_<time>.tsv` |
| `histograms` | `variables`: list; `bins`: default 10 | `FacetsDataset_histogram_<variable>_<time>.png` |
| `purity_cf` | | `purity_cf.txt` |
| `merged` | `files`: `adjusted_seg`, `unadjusted_seg`, `cncf`, `gene_level`, `qc`, `facets_qc`, `ccf_maf`, `nonsigned_maf` | `merged_*` |
| `cin` | `levels`: `full`, `chrom`, `arm`; `gene_or_seg`; `split_arms` | `cin_<level>.tsv` |
| `chromothripsis` | `gene_or_seg`, `min_breakpoints`, `min_ris_delta_pct`, `breakpoint_offset`, `bg_metric` | `chromothripsis_{cin_full,cin_arm,candidate_arms}.tsv` |
| `state_percents` | `include_x`, `include_y` | `hetloss_homdel_percent.tsv` |
| `ascets` | `genome`, `min_arm_breadth`, `arm_alt_frac_thresh`, `keep_noise`, `use_adjusted` | `ascets/` (needs R) |
| `copy_folders` | | `facets_folders/all/<P-00xxx>/<pair>/` |

Hypoploid and gain-heavy samples are flagged by `alteration`, in the `isHypoploidy` and `isGainHeavy` columns of `alteration_data.txt`: a sample is flagged when at least `alteration_arm_num` arms have at least `alteration_cov_min` percent altered.

`cin`, `chromothripsis`, `state_percents` and `ascets` need one fit per sample: `fits.single_run` true, and `source.fit_class` not `both`. `cin` and `chromothripsis` skip runs without a numeric ploidy and print how many.

Rendered scripts also write `run_summary.json` to the job folder, with run counts, build counts, unmatched sample IDs and the list of outputs.

### Repository safety

**Hard rule: a generated job never writes, modifies or deletes anything in a FACETS repository.**

- Methods that write into a repository (`backupMafs`, `FPTools`) and `persist_data` are not available from a spec.
- Every rendered script starts with `facetsAPI.write_guard`. For the rest of the Python process, writes are refused unless the real path, with symlinks resolved, is inside the job's output or job folder (or the temp folder) and outside every protected folder. Starting other programs is refused too, except the ASCETS R script.
- The ASCETS R code runs outside the guard. It has been read and writes only to its output folder. Any change to `tools/ascets` must keep that true.

### Command line

```bash
python -m facetsAPI.spec validate spec.json
python -m facetsAPI.spec render spec.json --output-dir results --job-dir . --facetsapi-home /path/to/facetsAPI > script.py
```

With a sample list, the script reads `samples.txt` from `--job-dir`. From Python, `spec.samples_text(spec)` returns its contents.

"""Chromothripsis candidate-arm flags from FACETS data.

Builds on the CIN scores in cin_scores.py: for every sample, each chromosome arm's CIN is compared with the
genome-wide CIN (relative instability score, RIS), its breakpoints are counted, and its cellular-fraction
states are checked for oscillation. An arm is flagged when ANY of:

  A) breakpoints >= min_breakpoints AND RIS >= background RIS * (1 + min_ris_delta_pct)
  B) breakpoints >= background breakpoints + breakpoint_offset
  C) oscillation: >= 2 CF states each occurring in >= 2 segments

One condition -> "suspect", two -> "highly suspect", three -> "critical suspect".

    python scripts/chromothripsis_scores.py --facets-dir <repo> --samples-file samples.txt -o out_prefix
writes <prefix>_cin_full.tsv, <prefix>_cin_arm.tsv and <prefix>_candidate_arms.tsv.
"""
import os
import sys
from collections import defaultdict, Counter
from statistics import mean, stdev, median

# facetsAPI must be importable: `pip install -e /path/to/facetsAPI` (editable) or add the repository
# root to PYTHONPATH. cin_scores.py is loaded from this scripts directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from facetsAPI import FacetsMeta, FacetsDataset  # noqa: E402
from cin_scores import cin_for_runs, isolated_run, write_cin_results  # noqa: E402


def identify_candidate_arms(
    cin_arm_results: list,
    full_genome_results: list,
    run_map: dict,
    min_breakpoints: int = 2,
    min_ris_delta_pct: float = 0.25,
    breakpoint_offset: int = 2,
    pseudocount: float = 1e-6,
    bg_metric: str = "median" #median or mean
) -> dict:
    """
    Flags candidate arms if ANY of:
      A) breakpoints ≥ min_breakpoints AND ris_raw ≥ avg_ris*(1+min_ris_delta_pct)
      B) breakpoints ≥ (avg_background_breaks_rounded + breakpoint_offset)
      D) oscillation (≥2 CF states each occurring ≥2 segments)

    RIS_raw = arm_CIN / genome_CIN  
    RIS_delta_pct = (RIS_raw / average_RIS_raw) - 1

    Returns per-sample:
      {
        "average_RIS_raw": float,
        "average_background_breakpoints": int,
        "candidates": [
          {
            "arm", "score", "breakpoints",
            "ris_raw", "ris_delta_pct",
            "met_conditions": [...],
            "cf_list", "n_unique_cf", "oscillation"
          }, …
        ]
      }
    """
    # map sample → genome-wide CIN
    full_cin = {r["sample"]: r["score"] for r in full_genome_results}

    all_results = {}

    for arm_res in cin_arm_results:
        sample     = arm_res["sample"]
        arms       = arm_res["arm"]
        arm_scores = arm_res["score_per_arm"]
        genome_cin = full_cin.get(sample, 0.0) or pseudocount

        # group segments by arm
        run = run_map[sample]
        run.defineArms()
        segs_by_arm = defaultdict(list)
        for seg in run.segments:
            segs_by_arm[seg.arm].append(seg)

        # compute ris_raw & breakpoints
        per_arm = []
        for arm, score in zip(arms, arm_scores):
            segs        = sorted(segs_by_arm.get(arm, []), key=lambda s: s.start)
            n_segs      = len(segs)
            breaks      = max(n_segs - 1, 0)
            ris_raw     = score / genome_cin

            per_arm.append({
                "arm":         arm,
                "score":       score,
                "breakpoints": breaks,
                "ris_raw":     ris_raw,
                "segments":    segs
            })

        # compute average RIS_raw
        ris_vals = [a["ris_raw"] for a in per_arm]
        if bg_metric == "median":
            avg_ris  = median(ris_vals) if ris_vals else 0.0
            mad = median(abs(x - median(ris_vals)) for x in ris_vals) or pseudocount

        if bg_metric == "mean":
            avg_ris  = mean(ris_vals) if ris_vals else 0.0

        #Compute stdev for z score calculation.
        if len(ris_vals) > 1:
            sd_ris = stdev(ris_vals)
        else:
            sd_ris = pseudocount

        # compute avg background breakpoints (rounded)
        bg_breaks     = [a["breakpoints"] for a in per_arm]
        if bg_metric == "median":
            avg_bg_breaks = round(median(bg_breaks)) if bg_breaks else 0
        if bg_metric == "mean":
            avg_bg_breaks = round(mean(bg_breaks)) if bg_breaks else 0

        # thresholds
        ris_thresh_A    = avg_ris * (1 + min_ris_delta_pct)
        breaks_thresh_B = avg_bg_breaks + breakpoint_offset

        # select candidates
        candidates = []
        for a in per_arm:
            # percent‐above‐average RIS
            ris_delta_pct = (a["ris_raw"] / avg_ris - 1.0) * 100  if avg_ris > 0 else 0.0

            if bg_metric == "median":
                z_ris = 0.6745 * (a["ris_raw"] - avg_ris) / mad
            if bg_metric == "mean":
                z_ris = (a["ris_raw"] - avg_ris) / (sd_ris or pseudocount)


            # build CF list & oscillation
            cf_list    = [round(seg.cf, 4) for seg in a["segments"]]
            counts     = Counter(cf_list)
            n_unique   = len(counts)
            osc_count  = sum(1 for cnt in counts.values() if cnt >= 2)
            oscillation = (osc_count >= 2)

            # test our three criteria
            conds = []
            if a["breakpoints"] >= min_breakpoints and a["ris_raw"] >= ris_thresh_A:
                conds.append("RIS Percent")
            if a["breakpoints"] >= breaks_thresh_B:
                conds.append("Breaks")
            if oscillation:
                conds.append("Oscillation")

            if not conds:
                continue

            # compute overall_flag
            n = len(conds)
            if n == 1:
                overall_flag = "suspect"
            elif n == 2:
                overall_flag = "highly suspect"
            else:  # n >= 3
                overall_flag = "critical suspect"

            candidates.append({
                "arm":             a["arm"],
                "score":           a["score"],
                "breakpoints":     a["breakpoints"],
                "ris_raw":         a["ris_raw"],
                "ris_delta_pct":   ris_delta_pct,
                "ris_zScore":      z_ris,
                "met_conditions":  conds,
                "overall_flag":    overall_flag,
                "cf_list":         cf_list,
                "n_unique_cf":     n_unique,
                "oscillation":     oscillation
            })


        all_results[sample] = {
            "average_RIS_raw":                 avg_ris,
            "average_background_breakpoints":  avg_bg_breaks,
            "candidates":                      candidates
        }

    return all_results

def arm_sort_key(arm_label: str):
    """
    Sort keys like '1p','1q','2p','2q',… so that numeric chromosome
    ascends and 'p' comes before 'q'.
    """
    label = str(arm_label)
    arm_letter = label[-1].lower()
    try:
        chrom = int(label[:-1])
    except ValueError:                      # X / Y / unlabelled arms sort after the autosomes
        chrom = {"x": 23, "y": 24}.get(label[:-1].lower(), 99)
    # p → 0, q → 1
    arm_index = 0 if arm_letter == "p" else 1
    return (chrom, arm_index)


def write_candidate_arms(summary_results: dict, output_file: str):
    """
    Write out a tab-delimited file with one row per candidate arm,
    sorted by sample then by natural arm order, now including z-scores.

    Columns:
      sample, arm, score, breakpoints,
      ris_raw, ris_delta_pct, ris_zScore,
      n_unique_cf, oscillation,
      RIS Percent, Breaks, Oscillation, overall_flag
    """
    # the three possible condition-names, updated:
    cond_columns = ["RIS Percent", "Breaks", "Oscillation"]

    with open(output_file, "w") as out:
        # write header
        headers = [
            "sample", "arm", "score", "breakpoints",
            "ris_raw", "ris_delta_pct", "ris_zScore",
            "n_unique_cf", "oscillation"
        ] + cond_columns + ["overall_flag"]
        out.write("\t".join(headers) + "\n")

        # iterate samples in sorted order
        for sample in sorted(summary_results):
            sample_info = summary_results[sample]
            candidates  = sample_info["candidates"]

            # sort arms naturally
            sorted_cands = sorted(candidates, key=lambda c: arm_sort_key(c["arm"]))

            for c in sorted_cands:
                # basic fields
                row = [
                    sample,
                    c["arm"],
                    f"{c['score']:.2f}",
                    str(c["breakpoints"]),
                    f"{c['ris_raw']:.4f}",
                    f"{c['ris_delta_pct']:.1f}",
                    f"{c.get('ris_zScore', 0):.2f}",
                    str(c["n_unique_cf"]),
                    str(c["oscillation"])
                ]
                # one column per condition
                for cond in cond_columns:
                    row.append("TRUE" if cond in c["met_conditions"] else "FALSE")

                # overall_flag last
                row.append(c["overall_flag"])

                out.write("\t".join(row) + "\n")


def chromothripsis_for_runs(runs, out_prefix, method="Eckert", gene_or_seg="seg", min_breakpoints=2,
                            min_ris_delta_pct=0.25, breakpoint_offset=2, bg_metric="median", verbose=False):
    """Score runs from an already built FacetsDataset, leaving them unchanged.

    Writes <out_prefix>_cin_full.tsv, <out_prefix>_cin_arm.tsv and <out_prefix>_candidate_arms.tsv.
    """
    common = dict(method=method, gene_or_seg=gene_or_seg, verbose=verbose)
    full_results = cin_for_runs(runs, splitArms=False, cinLevel="full", **common)
    arm_results = cin_for_runs(runs, splitArms=True, cinLevel="arm", **common)
    write_cin_results(full_results, out_prefix + "_cin_full.tsv")
    write_cin_results(arm_results, out_prefix + "_cin_arm.tsv")

    run_map = {run.id: isolated_run(run) for run in runs}
    summary = identify_candidate_arms(arm_results, full_results, run_map,
                                      min_breakpoints=min_breakpoints, min_ris_delta_pct=min_ris_delta_pct,
                                      breakpoint_offset=breakpoint_offset, bg_metric=bg_metric)
    write_candidate_arms(summary, out_prefix + "_candidate_arms.tsv")
    return summary


def main(argv=None):
    """CLI: genome-wide and per-arm CIN scores plus chromothripsis candidate-arm flags for a list of samples.

    Writes <out_prefix>_cin_full.tsv, <out_prefix>_cin_arm.tsv and <out_prefix>_candidate_arms.tsv.
    """
    import argparse
    ap = argparse.ArgumentParser(description="Chromothripsis candidate arms from FACETS: per-arm CIN relative to the genome, breakpoints, CF oscillation.")
    ap.add_argument("--facets-dir", required=True, help="FACETS repository root (the directory holding P-00xxx folders).")
    ap.add_argument("--clinical-sample-file", default="", help="data_clinical_sample(.oncokb).txt; optional.")
    ap.add_argument("--samples-file", help="One sample id per line (P-0001234-T01-IM6).")
    ap.add_argument("--sample", action="append", help="Repeatable; a sample id to process.")
    ap.add_argument("--method", default="Eckert", help="Label recorded in the output's method column.")
    ap.add_argument("--gene-or-seg", choices=("gene", "seg"), default="seg")
    ap.add_argument("--min-breakpoints", type=int, default=2)
    ap.add_argument("--min-ris-delta-pct", type=float, default=0.25)
    ap.add_argument("--breakpoint-offset", type=int, default=2)
    ap.add_argument("--bg-metric", choices=("median", "mean"), default="median")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output path prefix.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    samples = []
    if args.samples_file:
        with open(args.samples_file) as fh:
            samples.extend(ln.strip() for ln in fh if ln.strip() and not ln.startswith("#"))
    if args.sample:
        samples.extend(args.sample)
    if not samples:
        ap.error("no samples: use --samples-file or --sample")

    meta = FacetsMeta(args.clinical_sample_file, args.facets_dir, "purity")
    meta.setSingleRunPerSample(True, allowDefaults=True)
    meta.build_from_file_listing = True
    meta.samples_from_file.extend(samples)
    meta.buildFacetsMeta()
    ds = FacetsDataset(meta)
    ds.buildFacetsDataset()
    run_map = {run.id: run for run in ds.runList}
    present = [s for s in samples if s in run_map or any(k.startswith(s + "_") for k in run_map)]
    absent = [s for s in samples if s not in present]
    print("Samples with runs: %d; missing: %d" % (len(present), len(absent)))
    if absent:
        sys.stderr.write("[warn] no runs for: " + ", ".join(absent) + "\n")

    summary = chromothripsis_for_runs(ds.runList, args.out_prefix, method=args.method, gene_or_seg=args.gene_or_seg,
                                      min_breakpoints=args.min_breakpoints, min_ris_delta_pct=args.min_ris_delta_pct,
                                      breakpoint_offset=args.breakpoint_offset, bg_metric=args.bg_metric,
                                      verbose=args.verbose)
    print("Wrote %s_{cin_full,cin_arm,candidate_arms}.tsv" % args.out_prefix)
    return summary


if __name__ == "__main__":
    main()

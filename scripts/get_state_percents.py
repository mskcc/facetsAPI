#!/usr/bin/env python3
"""
Report genome-wide percent HETLOSS and HOMDEL per sample from FACETS segments.

Definitions (per request):
- HETLOSS: tcn == 1 and (lcn == 0 or lcn is None)   # lcn=-1 is treated as None
- HOMDEL : tcn == 0 and lcn == 0                     # strictly 0/0

Denominator:
- Total basepairs of all included chromosomes (autosomes by default).
- Use --include-x / --include-y to include sex chromosomes.

Output:
- TSV with one row per sample:
    sample  pct_HETLOSS  pct_HOMDEL  bp_HETLOSS  bp_HOMDEL  total_bp

Debug:
- --debug prints per-segment lines to stderr showing raw/norm tcn/lcn, bp,
  and whether the segment counted toward HETLOSS/HOMDEL/Neither.
"""

import sys
import os
import csv
import argparse
from collections import defaultdict

# facetsAPI must be importable: `pip install -e /path/to/facetsAPI` (editable) or add the repository
# root to PYTHONPATH.
from facetsAPI import FacetsMeta, FacetsDataset

# ---------- helpers ----------

def want_chrom(chrom, incX=False, incY=False):
    s = str(chrom).upper().replace('CHR', '')
    if s.isdigit():
        v = int(s)
        return 1 <= v <= 22
    if s == 'X':
        return incX
    if s == 'Y':
        return incY
    return False

def normalize_tcn_lcn(tcn, lcn):
    """
    Cast to ints when possible.
    Treat lcn < 0 (e.g., -1 sentinel) as None.
    If lcn > tcn, set lcn=None (invalid lesser copy).
    """
    try:
        tcn = int(round(float(tcn)))
    except Exception:
        tcn = None
    try:
        lcn = int(round(float(lcn)))
    except Exception:
        lcn = None

    if lcn is not None and lcn < 0:
        lcn = None
    if tcn is not None and lcn is not None and lcn > tcn:
        lcn = None
    return tcn, lcn

def classify_hetloss_homdel(tcn, lcn):
    """
    Returns one of: 'HETLOSS', 'HOMDEL', or None (for everything else),
    using the simplified rules described above.
    """
    tcnN, lcnN = normalize_tcn_lcn(tcn, lcn)

    # HOMDEL strictly 0/0
    if tcnN == 0 and lcnN == 0:
        return 'HOMDEL'

    # HETLOSS: 1/0 or 1/None
    if tcnN == 1 and (lcnN == 0 or lcnN is None):
        return 'HETLOSS'

    return None

def compute_sample_percentages(run, include_x=False, include_y=False, debug=False):
    """
    Walk all segments of this run and compute:
      total_bp, bp_HETLOSS, bp_HOMDEL, pct_HETLOSS, pct_HOMDEL
    """
    total_bp = 0
    bp_by_type = {'HETLOSS': 0, 'HOMDEL': 0}

    for seg in run.segments:
        if not want_chrom(seg.chrom, include_x, include_y):
            continue

        # segment length (inclusive)
        try:
            start = int(float(seg.start))
            end   = int(float(seg.end))
        except Exception:
            continue
        seg_len = max(0, end - start + 1)
        if seg_len == 0:
            continue

        # classification
        call = classify_hetloss_homdel(seg.tcn, seg.lcn)

        # update tallies
        total_bp += seg_len
        if call in bp_by_type:
            bp_by_type[call] += seg_len

        if debug:
            tcnN, lcnN = normalize_tcn_lcn(seg.tcn, seg.lcn)
            sys.stderr.write(
                f"[SEG] sample={run.id} chr={seg.chrom} {start}-{end} bp={seg_len} "
                f"raw_tcn={seg.tcn} raw_lcn={seg.lcn} norm_tcn={tcnN} norm_lcn={lcnN} -> {call or 'Neither'}\n"
            )

    pct_het = (bp_by_type['HETLOSS'] / total_bp * 100.0) if total_bp else 0.0
    pct_hom = (bp_by_type['HOMDEL']  / total_bp * 100.0) if total_bp else 0.0

    return {
        'total_bp': total_bp,
        'bp_HETLOSS': bp_by_type['HETLOSS'],
        'bp_HOMDEL':  bp_by_type['HOMDEL'],
        'pct_HETLOSS': pct_het,
        'pct_HOMDEL':  pct_hom,
    }

def write_summary(out_path, rows):
    """
    rows: list of dicts with keys:
      sample, pct_HETLOSS, pct_HOMDEL, bp_HETLOSS, bp_HOMDEL, total_bp
    """
    header = ["sample", "pct_HETLOSS", "pct_HOMDEL", "bp_HETLOSS", "bp_HOMDEL", "total_bp"]
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow([
                r["sample"],
                f'{r["pct_HETLOSS"]:.4f}',
                f'{r["pct_HOMDEL"]:.4f}',
                str(int(r["bp_HETLOSS"])),
                str(int(r["bp_HOMDEL"])),
                str(int(r["total_bp"])),
            ])

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Genome-wide percent HETLOSS and HOMDEL per sample from FACETS segments."
    )
    ap.add_argument("--facets-dir", required=True,
                    help="Base directory containing all FACETS output folders.")
    ap.add_argument("--clinical-sample-file", default="",
                    help="Path to clinical sample file (optional if using file listing).")
    ap.add_argument("--samples-file",
                    help="Text file with one sample/run ID per line (e.g., P-0001234-T01-IM6).")
    ap.add_argument("--sample", action="append",
                    help="Repeatable. Process this sample/run ID. (Alternative to --samples-file.)")
    ap.add_argument("--include-x", action="store_true", help="Include chromosome X in totals.")
    ap.add_argument("--include-y", action="store_true", help="Include chromosome Y in totals.")
    ap.add_argument("--debug", action="store_true",
                    help="Print per-segment classification details to stderr.")
    ap.add_argument("-o", "--output", required=True, help="Output TSV path.")
    args = ap.parse_args()

    # collect targets
    targets = []
    if args.samples_file:
        with open(args.samples_file) as f:
            targets.extend([ln.strip() for ln in f if ln.strip()])
    if args.sample:
        targets.extend(args.sample)
    if not targets:
        sys.exit("No samples provided. Use --samples-file or --sample.")

    # FACETS metadata/dataset
    meta = FacetsMeta(args.clinical_sample_file or "", args.facets_dir, "purity")
    meta.setSingleRunPerSample(True, allowDefaults=True)
    meta.build_from_file_listing = True
    meta.samples_from_file.extend(targets)
    meta.buildFacetsMeta()

    ds = FacetsDataset(meta)
    ds.buildFacetsDataset()

    # id -> run
    run_map = {run.id: run for run in ds.runList}

    out_rows = []
    missing = []

    for sample in targets:
        run = run_map.get(sample)
        if run is None:
            # allow prefix match like SAMPLE_* if your environment uses suffixed run IDs
            matches = [rid for rid in run_map.keys() if rid.startswith(sample + "_")]
            run = run_map.get(matches[0]) if matches else None

        if run is None:
            missing.append(sample)
            continue

        stats = compute_sample_percentages(run, include_x=args.include_x, include_y=args.include_y, debug=args.debug)
        out_rows.append({
            "sample": run.id,
            **stats
        })

    if missing:
        sys.stderr.write(f"[warn] Samples with no runs found: {', '.join(missing)}\n")

    write_summary(args.output, out_rows)

if __name__ == "__main__":
    main()

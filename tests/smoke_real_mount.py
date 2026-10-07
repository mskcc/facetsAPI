#!/usr/bin/env python3
"""Read-only smoke check of facetsAPI against a REAL FACETS repository (not run by unittest discover).

Builds a small sample list at runtime from the first few pair directories of one patient prefix,
then runs single-run and all-fits builds. Nothing is written to the repository; no ids are stored.

    python tests/smoke_real_mount.py /path/to/impact/facets/all/ [/path/to/impact_2n/facets/all/] [--prefix P-00000] [--n 3]

Exit status 0 when every repository given produced at least one run in both modes.
"""
import argparse
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from facetsAPI import FacetsMeta, FacetsDataset, MetaDictMap  # noqa: E402


def sample_ids(repo, prefix, n):
    pdir = os.path.join(repo, prefix)
    if not os.path.isdir(pdir):
        raise SystemExit("no such patient prefix dir: " + pdir)
    pairs = sorted(d for d in os.listdir(pdir) if os.path.isdir(os.path.join(pdir, d)))[:n]
    return [p.split("_")[0] for p in pairs]


def run(repo, prefix, n, fit_class):
    ids = sample_ids(repo, prefix, n)
    ok = True
    for single in (True, False):
        meta = FacetsMeta("", repo, "purity", fit_class=fit_class)
        meta.setSingleRunPerSample(single, single)
        meta.build_from_file_listing = True
        meta.samples_from_file.extend(ids)
        meta.setVerbose(True)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            meta.buildFacetsMeta()
            ds = FacetsDataset(meta)
            ds.buildFacetsDataset()
        mode = "single-run" if single else "all-fits"
        print("%s | %s | fit_class=%s | ids=%d runs=%d samples=%d" % (repo, mode, fit_class, len(ids), len(ds.runList), len(ds.sampleList)))
        print("   counts:", {k: v for k, v in meta.build_counts.items() if v})
        for key, entry in list(meta.master_file_dict.items())[:6]:
            print("   key shape:", "#" in key, "| class:", entry[MetaDictMap.FIT_CLASS], "| tier:", entry[MetaDictMap.FIT_TIER],
                  "| fit dir under repo:", entry[MetaDictMap.SELECTED_FIT_DIR].startswith(repo.rstrip("/")))
        if not ds.runList:
            ok = False
            print(buf.getvalue()[-3000:])
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("repos", nargs="+")
    ap.add_argument("--prefix", default="P-00000")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--fit-class", default="both")
    a = ap.parse_args()
    ok = all(run(r.rstrip("/") + "/", a.prefix, a.n, a.fit_class) for r in a.repos)
    print("SMOKE", "OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

"""Pre-generated (CADENCE) FacetsMeta: facetsAPI.prebuilt and the spec's fits.preset / source.prebuilt.

A prebuilt meta must be interchangeable with a live build of the same request: same runs, same fits,
same clinical data, same unmatched ids. Every test compares against the live build.
"""
import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from test_cadence_autoqc import build_repo  # noqa: E402

from facetsAPI import FacetsMeta, FitPolicy, MetaDictMap  # noqa: E402
from facetsAPI import spec as S  # noqa: E402
from facetsAPI.prebuilt import load_prebuilt_meta, find_prebuilt_meta  # noqa: E402

CHECKOUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def live_meta(repo, clinical, mode="purity", sample_file=None):
    m = FacetsMeta(clinical, repo, mode, fit_class="research", fit_policy=FitPolicy.cadence_autoqc(), keep_all_normals=True)
    m.setSingleRunPerSample(True, True)
    if sample_file:
        m.selectSamplesFromFile(sample_file)
    m.buildFacetsMeta()
    return m


def write_pickle(root, repo, clinical, stamp, mode="purity", **overrides):
    """Pickle a meta where BUILD_FACETS_META_2N would: <root>/cohort_level/<stamp>/research/autoqc/."""
    d = os.path.join(root, "cohort_level", stamp, "research", "autoqc")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "impact_facets_annotated.facets_meta.%s.pkl" % mode)
    m = FacetsMeta(clinical, repo, overrides.pop("built_mode", mode), path, fit_class="research",
                   fit_policy=overrides.pop("policy", FitPolicy.cadence_autoqc()),
                   keep_all_normals=overrides.pop("keep_all_normals", True))
    m.setSingleRunPerSample(True, True)
    m.cadence_built_at = stamp
    m.buildFacetsMeta()        # persist_data path does not exist -> builds and writes it
    return path


def view(meta):
    """What a job sees: run key -> (fit dir, fit status, cancer type)."""
    return {k: (v[MetaDictMap.SELECTED_FIT_DIR], v[MetaDictMap.FIT_STATUS], meta.cancer_type_map.get(meta.shortId(k)))
            for k, v in meta.master_file_dict.items()}


class PrebuiltLoaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.info = build_repo(self.root)
        self.repo, self.clin = self.info["repo"], self.info["clinical_file"]

    def tearDown(self):
        self.tmp.cleanup()

    def test_matches_live_build(self):
        write_pickle(self.root, self.repo, self.clin, "2026_11_10")
        meta, info = load_prebuilt_meta(self.repo, "purity", "research", self.clin)
        self.assertTrue(info["used"], info)
        self.assertEqual(info["stamp"], "2026_11_10")
        self.assertEqual(view(meta), view(live_meta(self.repo, self.clin)))
        self.assertEqual(meta.persist_data, "no")

    def test_sample_list_and_unmatched_match_live(self):
        write_pickle(self.root, self.repo, self.clin, "2026_11_10")
        tumor, pair1, _pair2 = self.info["multi"]
        other = self.info["pairs"]["autoqc_best"][1]
        ids = ["# header", tumor, other, "P-9999999-T01-IM6"]
        sf = os.path.join(self.root, "samples.txt")
        with open(sf, "w") as fh:
            fh.write("\n".join(ids) + "\n")
        meta, info = load_prebuilt_meta(self.repo, "purity", "research", self.clin, sf)
        live = live_meta(self.repo, self.clin, sample_file=sf)
        self.assertEqual(view(meta), view(live))
        self.assertEqual(meta.unmatched_sample_ids, live.unmatched_sample_ids)
        self.assertEqual(info["pairs_after_selection"], 3)

    def test_request_clinical_file_wins(self):
        write_pickle(self.root, self.repo, self.clin, "2026_11_10")
        tumor = self.info["multi"][0]
        dropped = self.info["pairs"]["autoqc_best"][0]
        lines = open(self.clin).read().splitlines(True)
        newer = os.path.join(self.root, "clin_newer.txt")
        with open(newer, "w") as fh:
            for line in lines:
                if line.startswith(dropped):
                    continue
                fh.write(line.replace("Lung Cancer", "Lung Cancer (new)") if line.startswith(tumor) else line)
        meta, _info = load_prebuilt_meta(self.repo, "purity", "research", newer)
        self.assertEqual(view(meta), view(live_meta(self.repo, newer)))
        self.assertNotIn(dropped, {meta.shortId(k) for k in meta.master_file_dict})
        self.assertEqual(meta.cancer_type_map[tumor], "Lung Cancer (new)")
        self.assertEqual(meta.clinical_sample_file, newer)

    def test_newest_usable_stamp_wins_and_mismatches_are_skipped(self):
        write_pickle(self.root, self.repo, self.clin, "2026_11_03")
        write_pickle(self.root, self.repo, self.clin, "2026_11_10")
        # newest is built with the stock policy and only the last normal: must not be used
        write_pickle(self.root, self.repo, self.clin, "2026_11_17", policy=FitPolicy(), keep_all_normals=False)
        os.makedirs(os.path.join(self.root, "cohort_level", "notes"))
        self.assertEqual([s for s, _p in find_prebuilt_meta(self.repo, "purity")], ["2026_11_17", "2026_11_10", "2026_11_03"])
        meta, info = load_prebuilt_meta(self.repo, "purity", "research", self.clin)
        self.assertEqual(info["stamp"], "2026_11_10")
        self.assertIn("2026_11_17: fit_policy", info["skipped_newer"][0])

    def test_falls_back_when_nothing_fits(self):
        meta, info = load_prebuilt_meta(self.repo, "purity", "research", self.clin)
        self.assertIsNone(meta)
        self.assertIn("no purity pickle", info["reason"])
        write_pickle(self.root, self.repo, self.clin, "2026_11_10")
        self.assertIsNone(load_prebuilt_meta(self.repo, "hisens", "research", self.clin)[0])
        meta, info = load_prebuilt_meta(self.repo, "purity", "clinical", self.clin)
        self.assertIsNone(meta)
        self.assertIn("only for fit_class=research", info["reason"])
        # a pickle of another repository is not this one's
        other = os.path.join(self.root, "elsewhere", "all")
        os.makedirs(other)
        shutil.copytree(os.path.join(self.root, "cohort_level"), os.path.join(self.root, "elsewhere", "cohort_level"))
        meta, info = load_prebuilt_meta(other, "purity", "research", self.clin)
        self.assertIsNone(meta)
        self.assertIn("facets_repo_path", info["reason"])

    def test_wrong_mode_inside_file_and_garbage_are_skipped(self):
        write_pickle(self.root, self.repo, self.clin, "2026_11_10", mode="hisens", built_mode="purity")
        d = os.path.join(self.root, "cohort_level", "2026_11_17", "research", "autoqc")
        os.makedirs(d)
        with open(os.path.join(d, "x.facets_meta.hisens.pkl"), "wb") as fh:
            fh.write(b"not a pickle")
        meta, info = load_prebuilt_meta(self.repo, "hisens", "research", self.clin)
        self.assertIsNone(meta)
        self.assertIn("2026_11_17: unreadable", info["reason"])
        self.assertIn("2026_11_10: hisens_vs_purity", info["reason"])


def spec(**over):
    s = {"spec_version": S.SPEC_VERSION, "name": "pre", "samples": {"source": "all"}, "filters": {},
         "outputs": [{"type": "report"}],
         "source": {"facets_dir": "/repo/all/", "clinical_file": "/repo/clin.txt", "mode": "purity",
                    "fit_class": "research", "prebuilt": True},
         "fits": {"single_run": True, "preset": "cadence_autoqc"}}
    for k, v in over.items():
        s[k] = v
    return s


class SpecTests(unittest.TestCase):
    def test_validation(self):
        self.assertEqual(S.validate(spec())[0], [])
        bad = copy.deepcopy(spec())
        bad["fits"] = {"single_run": True, "tiers": ["autoqc", "default"]}
        self.assertTrue(any("Pre-generated" in e for e in S.validate(bad)[0]))
        bad = copy.deepcopy(spec())
        bad["source"]["fit_class"] = "clinical"
        self.assertTrue(any("Pre-generated" in e for e in S.validate(bad)[0]))
        bad = copy.deepcopy(spec())
        bad["source"]["prebuilt"] = "yes"
        self.assertIn("source.prebuilt must be true or false.", S.validate(bad)[0])
        bad = copy.deepcopy(spec())
        bad["fits"] = {"single_run": False, "preset": "cadence_autoqc"}
        self.assertIn("fits.preset needs fits.single_run true.", S.validate(bad)[0])
        bad["fits"] = {"single_run": True, "preset": "nope"}
        self.assertTrue(any(e.startswith("fits.preset must be one of") for e in S.validate(bad)[0]))
        ok = copy.deepcopy(spec())
        ok["source"]["prebuilt"] = False
        ok["source"]["fit_class"] = "both"
        self.assertEqual(S.validate(ok)[0], [])      # the preset alone works for any class

    def test_render_compiles(self):
        c = {"job_id": "j", "user": "u", "facetsapi_home": CHECKOUT, "facetsapi_sha": "x", "output_dir": "/o", "job_dir": "/j"}
        for prebuilt in (True, False):
            s = copy.deepcopy(spec())
            s["source"]["prebuilt"] = prebuilt
            code = S.render(s, c)
            compile(code, "script.py", "exec")
            self.assertIn("FitPolicy.cadence_autoqc(), keep_all_normals=True", code)
            self.assertEqual("load_prebuilt_meta" in code, prebuilt)


class RenderedScriptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.info = build_repo(self.root)
        self.n = 0

    def tearDown(self):
        self.tmp.cleanup()

    def snapshot(self):
        snap = {}
        for base in (self.info["repo"], os.path.join(self.root, "cohort_level")):
            for r, ds, fs in os.walk(base):
                for name in ds + fs:
                    p = os.path.join(r, name)
                    st = os.lstat(p)
                    digest = hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.isfile(p) else None
                    snap[p] = (st.st_mode, st.st_size, st.st_mtime_ns, digest)
        return snap

    def run_spec(self, s, samples=None):
        self.n += 1
        job = os.path.join(self.root, "job%d" % self.n)
        out = os.path.join(job, "out")
        os.makedirs(job)
        s = copy.deepcopy(s)
        s["source"]["facets_dir"] = self.info["repo"]
        s["source"]["clinical_file"] = self.info["clinical_file"]
        if samples:
            s["samples"] = {"source": "list", "ids": samples}
            with open(os.path.join(job, "samples.txt"), "w") as fh:
                fh.write(S.samples_text(s))
        script = os.path.join(job, "script.py")
        protected = [self.info["repo"], os.path.join(self.root, "cohort_level")]
        with open(script, "w") as fh:
            fh.write(S.render(s, {"job_id": "j", "user": "u", "facetsapi_home": CHECKOUT, "facetsapi_sha": "x",
                                  "output_dir": out, "job_dir": job, "protected_roots": protected}))
        r = subprocess.run([sys.executable, script], capture_output=True, text=True, cwd=job)
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        with open(os.path.join(job, "run_summary.json")) as fh:
            return json.load(fh), open(os.path.join(out, "report.txt")).read(), r.stdout

    def test_prebuilt_and_live_jobs_agree_and_nothing_is_written(self):
        write_pickle(self.root, self.info["repo"], self.info["clinical_file"], "2026_11_10")
        before = self.snapshot()
        pre, pre_report, stdout = self.run_spec(spec())
        live_spec = spec()
        live_spec["source"]["prebuilt"] = False
        live, live_report, _ = self.run_spec(live_spec)
        self.assertEqual(self.snapshot(), before)
        self.assertTrue(pre["meta_source"]["used"])
        self.assertIn("Using pre-generated FacetsMeta", stdout)
        self.assertEqual(live["meta_source"], {"requested": False, "used": False})
        self.assertEqual(pre["num_runs"], live["num_runs"])
        self.assertEqual(pre_report, live_report)

    def test_sample_list(self):
        write_pickle(self.root, self.info["repo"], self.info["clinical_file"], "2026_11_10")
        ids = [self.info["multi"][0], "P-9999999-T01-IM6"]
        pre, pre_report, _ = self.run_spec(spec(), samples=ids)
        live_spec = spec()
        live_spec["source"]["prebuilt"] = False
        live, live_report, _ = self.run_spec(live_spec, samples=ids)
        self.assertEqual((pre["num_runs"], pre["unmatched_sample_ids"]), (live["num_runs"], live["unmatched_sample_ids"]))
        self.assertEqual(pre["num_runs"], 2)
        self.assertEqual(pre_report, live_report)

    def test_falls_back_to_a_live_build(self):
        summary, _report, stdout = self.run_spec(spec())
        self.assertFalse(summary["meta_source"]["used"])
        self.assertIn("no purity pickle", summary["meta_source"]["reason"])
        self.assertIn("Building it from the repository", stdout)
        self.assertGreater(summary["num_runs"], 0)


if __name__ == "__main__":
    unittest.main()

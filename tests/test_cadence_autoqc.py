"""FitPolicy.cadence_autoqc() and FacetsMeta(keep_all_normals=True).

cadence_autoqc() must reproduce CADENCE's 2n cohort policy 'autoqc'
(facets-preview compile_cohort_annotations_2n): the 'auto-qc' reviewer's
auto_qc_best_fit, else the default fit -- ignoring human verdicts, auto_qc_pass,
legacy auto-qc-script rows and ultra fits, with no fall-through. Each fixture
pair below is one way the stock policy and CADENCE disagree.
"""
import glob
import os
import pickle
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import write_sample_dir, _pair, _w  # noqa: E402

from facetsAPI import FacetsMeta, FacetsDataset, MetaDictMap  # noqa: E402
from facetsAPI.exceptions import FacetsConfigError  # noqa: E402
from facetsAPI.fitpolicy import FitPolicy, classify  # noqa: E402
from facetsAPI.layout import read_manifest  # noqa: E402

AQ = "auto-qc"

# name -> (pair number, research fits, CADENCE-expected fit, stock-policy fit)
CASES = {
    # no best fit; a non-default fit passed -> CADENCE default, stock takes the pass
    "pass_only": (101, [
        {"fit_name": "default", "status": "not_reviewed"},
        {"fit_name": "alt_dipLogR_0.2", "status": "auto_qc_pass", "reviewed_by": AQ},
    ], "default", "alt_dipLogR_0.2"),
    # human best fit + autoQC best on an ultra fit -> CADENCE default; stock takes the human fit
    "human_and_ultra_best": (102, [
        {"fit_name": "default", "status": "auto_qc_pass", "reviewed_by": AQ},
        {"fit_name": "alt_dipLogR_-0.1", "status": "reviewed_best_fit", "reviewed_by": "alice"},
        {"fit_name": "alt_dipLogR_0.3/ultra", "status": "auto_qc_best_fit", "reviewed_by": AQ, "ultra": True},
    ], "default", "alt_dipLogR_-0.1"),
    # ordinary autoQC best fit -> both agree
    "autoqc_best": (103, [
        {"fit_name": "default", "status": "auto_qc_pass", "reviewed_by": AQ},
        {"fit_name": "alt_dipLogR_0.1", "status": "auto_qc_best_fit", "reviewed_by": AQ},
    ], "alt_dipLogR_0.1", "alt_dipLogR_0.1"),
    # legacy autoQC (reviewed_best_fit by *auto-qc-script) -> CADENCE ignores it
    "legacy_autoqc": (104, [
        {"fit_name": "default", "status": "not_reviewed"},
        {"fit_name": "alt_dipLogR_0.5", "status": "reviewed_best_fit", "reviewed_by": "richara4/auto-qc-script"},
    ], "default", "alt_dipLogR_0.5"),
    # auto_qc_best_fit written by a different reviewer -> CADENCE requires reviewed_by == 'auto-qc'
    "other_reviewer": (105, [
        {"fit_name": "default", "status": "not_reviewed"},
        {"fit_name": "alt_dipLogR_0.4", "status": "auto_qc_best_fit", "reviewed_by": "someone-else"},
    ], "default", "alt_dipLogR_0.4"),
    # best fit lacks purity files -> CADENCE does not fall through (skip); stock falls to default
    "best_missing_files": (106, [
        {"fit_name": "default", "status": "not_reviewed"},
        {"fit_name": "alt_dipLogR_0.7", "status": "auto_qc_best_fit", "reviewed_by": AQ},
    ], None, "default"),
}


def build_repo(root):
    repo = os.path.join(root, "all")
    info = {"repo": repo + "/", "pairs": {}}

    def write_2n(tumor, pair, research_fits):
        base = os.path.join(repo, tumor[0:7], pair)
        write_sample_dir(os.path.join(base, "clinical"), pair, [{"fit_name": "default", "status": "not_reviewed"}])
        write_sample_dir(os.path.join(base, "research"), pair, research_fits)
        return base

    for name, (n, fits, _cad, _stock) in CASES.items():
        tumor, pair = _pair(n)
        base = write_2n(tumor, pair, fits)
        info["pairs"][name] = (tumor, pair)
        if name == "best_missing_files":
            for f in glob.glob(os.path.join(base, "research", "alt_dipLogR_0.7", "*_purity*")):
                os.remove(f)

    # one tumor, two normals, both 2n: N01 has an autoQC best alt fit, N02 only default
    tumor, pair1 = _pair(107)
    _, pair2 = _pair(107, nn="N02")
    write_2n(tumor, pair1, [{"fit_name": "default", "status": "not_reviewed"},
                            {"fit_name": "alt_dipLogR_0.2", "status": "auto_qc_best_fit", "reviewed_by": AQ}])
    write_2n(tumor, pair2, [{"fit_name": "default", "status": "not_reviewed"}])
    info["multi"] = (tumor, pair1, pair2)

    tumors = sorted({t for t, _p in info["pairs"].values()} | {tumor})
    rows = ["SAMPLE_ID\tPATIENT_ID\tCANCER_TYPE\tCANCER_TYPE_DETAILED\tTUMOR_PURITY\tONCOTREE_CODE\tCVR_TMB_SCORE\tMSI_SCORE\n"]
    for i, t in enumerate(tumors):
        ct = "Lung Cancer" if t == tumor else "Breast Cancer"
        rows.append("\t".join([t, t[0:9], ct, ct + " NOS", "50", "LUAD" if t == tumor else "IDC", "3.0", "0.5"]) + "\n")
    info["clinical_file"] = os.path.join(root, "data_clinical_sample.oncokb.txt")
    _w(info["clinical_file"], "".join(rows))
    return info


class CadencePolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build_repo(cls.tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def research(self, name):
        tumor, pair = self.info["pairs"][name]
        d = os.path.join(self.info["repo"], tumor[0:7], pair, "research")
        return read_manifest(os.path.join(d, "facets_review.manifest")), d

    def test_select_fit_matches_cadence_and_documents_stock_difference(self):
        for name, (_n, _fits, cadence_fit, stock_fit) in CASES.items():
            df, d = self.research(name)
            with self.subTest(case=name):
                cand, _files, _tried = FitPolicy.cadence_autoqc().select_fit(df, "purity", d)
                self.assertEqual(cand.fit_name if cand else None, cadence_fit)
                stock, _f, _t = FitPolicy(("autoqc", "default")).select_fit(df, "purity", d)
                if name in ("human_and_ultra_best",):
                    stock, _f, _t = FitPolicy().select_fit(df, "purity", d)  # reviewed > autoqc > default
                self.assertEqual(stock.fit_name if stock else None, stock_fit)

    def test_labels_and_tiers(self):
        df, d = self.research("autoqc_best")
        cand, _f, _t = FitPolicy.cadence_autoqc().select_fit(df, "purity", d)
        self.assertEqual((cand.tier, cand.label), ("autoqc", "AutoQC_Best"))
        self.assertEqual(cand.fit_status, "auto_qc_best_fit#alt_dipLogR_0.1#auto-qc")
        df, d = self.research("pass_only")
        cand, _f, _t = FitPolicy.cadence_autoqc().select_fit(df, "purity", d)
        self.assertEqual((cand.fit_name, cand.tier, cand.label), ("default", "default", "Default"))
        # excluded rows keep their status label; a default-named one lands in the default tier
        self.assertEqual(classify("auto_qc_pass", "default", AQ, autoqc_statuses=("auto_qc_best_fit",)),
                         ("default", 0, "AutoQC_Pass"))
        self.assertEqual(classify("auto_qc_pass", "alt_dipLogR_0.2", AQ, autoqc_statuses=("auto_qc_best_fit",)),
                         ("unreviewed", 0, "AutoQC_Pass"))
        self.assertEqual(classify("reviewed_best_fit", "alt", "x/auto-qc-script", legacy_autoqc=False),
                         ("unreviewed", 0, "AutoQC_Best"))
        self.assertEqual(classify("auto_qc_best_fit", "alt", "other", autoqc_reviewer=AQ), ("unreviewed", 0, "AutoQC_Best"))
        # defaults unchanged
        self.assertEqual(classify("auto_qc_pass", "alt", AQ), ("autoqc", 1, "AutoQC_Pass"))
        self.assertEqual(classify("reviewed_best_fit", "alt", "x/auto-qc-script"), ("autoqc", 0, "AutoQC_Best"))

    def test_validation_and_with_defaults(self):
        with self.assertRaises(FacetsConfigError):
            FitPolicy(autoqc_statuses=())
        with self.assertRaises(FacetsConfigError):
            FitPolicy(autoqc_statuses=("reviewed_best_fit",))
        with self.assertRaises(FacetsConfigError):
            FitPolicy(autoqc_statuses=("auto_qc_pass", "auto_qc_pass"))
        with self.assertRaises(FacetsConfigError):
            FitPolicy(autoqc_reviewer="")
        self.assertEqual(FitPolicy(autoqc_statuses="auto_qc_pass").autoqc_statuses, ("auto_qc_pass",))
        p = FitPolicy.cadence_autoqc().with_defaults(False)
        self.assertEqual(p.tiers, ("autoqc",))
        self.assertEqual((p.autoqc_statuses, p.legacy_autoqc, p.autoqc_reviewer, p.fall_through),
                         (("auto_qc_best_fit",), False, AQ, False))
        self.assertEqual(FitPolicy.cadence_autoqc(include_ultra_fits=True).include_ultra_fits, True)

    def test_old_pickled_policy_gets_new_field_defaults(self):
        old = FitPolicy()
        state = dict(old.__dict__)
        for k in ("autoqc_statuses", "legacy_autoqc", "autoqc_reviewer", "fall_through"):
            state.pop(k)
        revived = FitPolicy.__new__(FitPolicy)
        revived.__setstate__(state)
        self.assertEqual(revived, FitPolicy())
        self.assertEqual(pickle.loads(pickle.dumps(FitPolicy.cadence_autoqc())), FitPolicy.cadence_autoqc())


class KeepAllNormalsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build_repo(cls.tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def meta(self, keep, persist="no", mode="purity"):
        m = FacetsMeta(self.info["clinical_file"], self.info["repo"], mode, persist, fit_class="research",
                       fit_policy=FitPolicy.cadence_autoqc(), keep_all_normals=keep)
        m.setSingleRunPerSample(True, True)
        m.buildFacetsMeta()
        return m

    def test_pair_keys_and_cadence_selection(self):
        m = self.meta(True)
        tumor, pair1, pair2 = self.info["multi"]
        expected = {self.info["pairs"][n][1]: f for n, (_x, _y, f, _s) in CASES.items() if f is not None}
        expected[pair1] = "alt_dipLogR_0.2"
        expected[pair2] = "default"
        self.assertEqual(set(m.master_file_dict), set(expected))
        for key, fit in expected.items():
            self.assertTrue(m.master_file_dict[key][MetaDictMap.SELECTED_FIT_DIR].endswith("/research/%s/" % fit), key)
            self.assertEqual(m.pairIdForKey(key), key)
            self.assertEqual(m.shortId(key), key.split("_")[0])
            self.assertEqual(m.master_file_dict[key][MetaDictMap.FIT_CLASS], "research")
        self.assertEqual(m.long_id_map[pair2], [pair2])
        self.assertEqual(m.build_counts["multi_normal"], 1)
        self.assertEqual(m.build_counts["missing_files"], 1)        # best_missing_files skipped, not defaulted
        self.assertEqual(m.build_counts["AutoQC_Best"], 2)
        # facetsAPI convention: the label is the row's own status, the tier is where the policy found it --
        # human_and_ultra_best's default row is an auto_qc_pass row taken as the DEFAULT fallback.
        hu = self.info["pairs"]["human_and_ultra_best"][1]
        self.assertEqual(m.fit_map[hu][0], "AutoQC_Pass")
        self.assertEqual(m.master_file_dict[hu][MetaDictMap.FIT_TIER], "default")
        self.assertEqual(m.build_counts["Default"] + m.build_counts["AutoQC_Pass"], len(expected) - 2)
        tiers = [v[MetaDictMap.FIT_TIER] for v in m.master_file_dict.values()]
        self.assertEqual((tiers.count("autoqc"), tiers.count("default")), (2, len(expected) - 2))

    def test_default_keeps_last_normal_only(self):
        m = self.meta(False)
        tumor, pair1, pair2 = self.info["multi"]
        self.assertIn(tumor, m.master_file_dict)
        self.assertNotIn(pair1, m.master_file_dict)
        self.assertTrue(m.master_file_dict[tumor][MetaDictMap.SELECTED_FIT_DIR].endswith(pair2 + "/research/default/"))
        self.assertEqual(m.pairIdForKey(tumor), pair2)

    def test_dataset_builds_with_clinical_filters_on_pair_keys(self):
        m = self.meta(True)
        tumor, pair1, pair2 = self.info["multi"]
        ds = FacetsDataset(m)
        ds.setCancerTypeFilter(["Lung Cancer"])
        ds.buildFacetsDataset()
        self.assertEqual(sorted(r.id for r in ds.runList), sorted([pair1, pair2]))
        self.assertEqual({r.cancerType for r in ds.runList}, {"Lung Cancer"})

    def test_persist_round_trip(self):
        path = os.path.join(self.tmp.name, "meta.pkl")
        built = self.meta(True, persist=path)
        self.assertTrue(os.path.isfile(path))
        loaded = FacetsMeta("", self.info["repo"], "purity", path)
        loaded.buildFacetsMeta()
        self.assertEqual(set(loaded.master_file_dict), set(built.master_file_dict))
        self.assertEqual(loaded.fit_policy, FitPolicy.cadence_autoqc())
        self.assertTrue(loaded.keep_all_normals and loaded.selectSingleRun)
        os.remove(path)

    def test_hisens_mode(self):
        m = self.meta(True, mode="hisens")
        _t, pair1, _p2 = self.info["multi"]
        self.assertTrue(m.master_file_dict[pair1][MetaDictMap.OUT_FILE].endswith("_hisens.out"))
        # hisens files exist for best_missing_files' alt fit, so it is selected there
        bm = self.info["pairs"]["best_missing_files"][1]
        self.assertTrue(m.master_file_dict[bm][MetaDictMap.SELECTED_FIT_DIR].endswith("/research/alt_dipLogR_0.7/"))


if __name__ == "__main__":
    unittest.main()

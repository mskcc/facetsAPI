import contextlib
import io
import os
import pickle
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

from facetsAPI import FacetsMeta, FacetsDataset, MetaDictMap  # noqa: E402
from facetsAPI.exceptions import FacetsConfigError, FacetsRepositoryError  # noqa: E402
from facetsAPI.fitpolicy import FitPolicy  # noqa: E402


class _OldPickledMeta:
    """Stand-in for a FacetsMeta pickled by facetsAPI 0.5 (no _persist_version attribute)."""


def quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


class MetaDatasetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.repo = cls.info["repo"]
        cls.clin = cls.info["clinical_file"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def meta(self, single=True, defaults=True, fit_class="research", clinical=True, policy=None, mode="purity"):
        m = FacetsMeta(self.clin if clinical else "", self.repo, mode, fit_class=fit_class, fit_policy=policy)
        m.setSingleRunPerSample(single, defaults if single else False)
        quiet(m.buildFacetsMeta)
        return m

    # ---- construction / validation ----------------------------------------------------
    def test_constructor_validation(self):
        with self.assertRaises(FacetsConfigError):
            FacetsMeta(self.clin, self.repo, "purity", fit_class="ultra")
        with self.assertRaises(FacetsConfigError):
            FacetsMeta(self.clin, self.repo, "hisense")
        with self.assertRaises(FacetsConfigError):
            FacetsMeta(self.clin, self.repo, "purity", fit_policy="reviewed")
        m = FacetsMeta(self.clin, "", "purity")
        with self.assertRaises(FacetsRepositoryError):
            quiet(m.buildFacetsMeta)
        m = FacetsMeta(self.clin, os.path.join(self.tmp.name, "nope"), "purity")
        with self.assertRaises(FacetsRepositoryError):
            quiet(m.buildFacetsMeta)
        m = FacetsMeta(os.path.join(self.tmp.name, "missing.txt"), self.repo, "purity")
        with self.assertRaises(FacetsRepositoryError):
            quiet(m.buildFacetsMeta)
        with self.assertRaises(FacetsConfigError):
            FacetsMeta(self.clin, self.repo).setSingleRunPerSample(False, True)
        with self.assertRaises(FacetsConfigError):
            FacetsMeta(self.clin, self.repo).setSingleRunPerSample("yes", True)

    def test_state_is_per_instance(self):
        a = FacetsMeta(self.clin, self.repo, "hisens")
        b = FacetsMeta(self.clin, self.repo, "purity")
        a.setSingleRunPerSample(True, False)
        self.assertTrue(a.selectSingleRun)
        self.assertFalse(b.selectSingleRun)
        self.assertEqual(a.hisens_vs_purity, "hisens")
        self.assertEqual(b.hisens_vs_purity, "purity")
        a.failed_samples.append("x")
        self.assertEqual(b.failed_samples, [])
        self.assertEqual(FacetsMeta.failed_samples, [])     # class default untouched
        a.setVerbose(True)
        self.assertTrue(a.run_verbose)
        self.assertFalse(b.run_verbose)

    def test_policy_mapping_from_setSingleRunPerSample(self):
        m = FacetsMeta(self.clin, self.repo)
        m.setSingleRunPerSample(True, True)
        self.assertEqual(m.fit_policy.tiers, ("reviewed", "autoqc", "default"))
        m.setSingleRunPerSample(True, False)
        self.assertEqual(m.fit_policy.tiers, ("reviewed", "autoqc"))
        # an explicit policy survives setSingleRunPerSample, only the default tier toggles
        m = FacetsMeta(self.clin, self.repo, fit_policy=FitPolicy.autoqc_only())
        m.setSingleRunPerSample(True, True)
        self.assertEqual(m.fit_policy.tiers, ("autoqc", "default"))
        m.setSingleRunPerSample(True, False)
        self.assertEqual(m.fit_policy.tiers, ("autoqc",))
        m.setFitPolicy(FitPolicy.reviewed_only())
        self.assertEqual(m.fit_policy.tiers, ("reviewed",))
        self.assertFalse(m.use_unreviewed_defaults)

    # ---- single-run selection ----------------------------------------------------------
    def test_single_run_default_policy_standard(self):
        m = self.meta()
        std = self.info["standard"]
        keys = m.master_file_dict
        self.assertEqual(keys[std["reviewed_best"][0]][MetaDictMap.FIT_STATUS], "reviewed_best_fit#alt_dipLogR_-0.1#alice")
        self.assertEqual(keys[std["autoqc_only"][0]][MetaDictMap.FIT_STATUS], "auto_qc_best_fit#alt_dipLogR_0.2#auto-qc")
        self.assertEqual(keys[std["default_only"][0]][MetaDictMap.FIT_STATUS], "not_reviewed#default#NA")
        self.assertEqual(keys[std["no_comment_acceptable"][0]][MetaDictMap.FIT_STATUS], "reviewed_acceptable_fit#default#carol")
        self.assertNotIn(std["no_manifest"][0], keys)
        self.assertNotIn(std["no_dir"][0], keys)
        self.assertEqual(keys[std["reviewed_best"][0]][MetaDictMap.FIT_TIER], "reviewed")
        self.assertEqual(keys[std["autoqc_only"][0]][MetaDictMap.FIT_TIER], "autoqc")
        self.assertEqual(keys[std["default_only"][0]][MetaDictMap.FIT_TIER], "default")
        self.assertEqual(keys[std["reviewed_best"][0]][MetaDictMap.FIT_CLASS], "standard")
        self.assertEqual(m.fit_map[std["reviewed_best"][0]][0], "Best")
        self.assertEqual(m.fit_map[std["autoqc_only"][0]][0], "AutoQC_Best")
        self.assertEqual(m.fit_map[std["default_only"][0]][0], "Default")
        c = m.build_counts
        self.assertEqual((c["Best"], c["Acceptable"], c["AutoQC_Best"], c["Default"]), (2, 1, 3, 2))   # P-0000010: reviewed_best_fit by an auto-qc-script reviewer counts as autoQC
        self.assertEqual((c["not_found"], c["missing_manifest"], c["multi_normal"]), (1, 1, 1))

    def test_fit_dirs_resolve_from_discovered_dir_not_manifest_path(self):
        m = self.meta()
        for key, entry in m.master_file_dict.items():
            fit_dir = entry[MetaDictMap.SELECTED_FIT_DIR]
            self.assertTrue(fit_dir.startswith(self.repo), fit_dir)
            self.assertFalse(fit_dir.startswith("/nonexistent"), fit_dir)
            self.assertTrue(fit_dir.endswith("/"))
            self.assertTrue(os.path.isdir(fit_dir))
            for slot in (MetaDictMap.OUT_FILE, MetaDictMap.CNCF_FILE, MetaDictMap.QC_FILE, MetaDictMap.GENE_FILE,
                         MetaDictMap.ADJUSTED_SEG_FILE, MetaDictMap.UNADJUSTED_SEG_FILE, MetaDictMap.FACETS_QC_FILE,
                         MetaDictMap.MANIFEST_FILE):
                self.assertTrue(os.path.isfile(entry[slot]), (key, slot, entry[slot]))
            self.assertTrue(os.path.isdir(entry[MetaDictMap.SAMPLE_BASE_DIR]))
        # Rdata optional
        self.assertEqual(m.master_file_dict[self.info["standard"]["no_comment_acceptable"][0]][MetaDictMap.RDATA_FILE], "")

    def test_multi_normal_deterministic(self):
        m = self.meta()
        tumor, pair1, pair2 = self.info["standard"]["multi_normal"]
        self.assertIn(pair2, m.master_file_dict[tumor][MetaDictMap.SELECTED_FIT_DIR])
        self.assertEqual(m.long_id_map[tumor], [pair2])

    def test_only_policies_on_meta(self):
        std = self.info["standard"]
        m = self.meta(policy=FitPolicy.reviewed_only(), defaults=False)
        self.assertIn(std["reviewed_best"][0], m.master_file_dict)
        self.assertIn(std["no_comment_acceptable"][0], m.master_file_dict)
        self.assertNotIn(std["autoqc_only"][0], m.master_file_dict)
        self.assertNotIn(std["default_only"][0], m.master_file_dict)
        self.assertEqual(m.build_counts["no_fit_selected"], 5)   # autoqc-only, default-only, the auto-reviewed one, and both 2n research classes
        m = self.meta(policy=FitPolicy.autoqc_only(), defaults=False)
        self.assertEqual(sorted(k[0:9] for k in m.master_file_dict), ["P-0000002", "P-0000007", "P-0000010"])
        m = self.meta(policy=FitPolicy.default_only(), defaults=True)
        self.assertNotIn(std["reviewed_best"][0] + "x", m.master_file_dict)
        self.assertEqual(m.master_file_dict[std["reviewed_best"][0]][MetaDictMap.FIT_STATUS], "not_reviewed#default#NA")
        self.assertEqual(m.master_file_dict[std["autoqc_only"][0]][MetaDictMap.FIT_STATUS], "auto_qc_pass#default#auto-qc")
        # legacy flag: reviewed+autoqc without defaults
        m = self.meta(defaults=False)
        self.assertNotIn(std["default_only"][0], m.master_file_dict)
        self.assertIn(std["autoqc_only"][0], m.master_file_dict)

    # ---- 2n ----------------------------------------------------------------------------
    def test_2n_fit_class_research_default(self):
        m = self.meta()
        tumor, pair = self.info["twon"]["full"]
        e = m.master_file_dict[tumor]
        self.assertEqual(e[MetaDictMap.FIT_CLASS], "research")
        self.assertTrue(e[MetaDictMap.SAMPLE_BASE_DIR].endswith(pair + "/research"))
        self.assertTrue(e[MetaDictMap.FACETS_QC_FILE].endswith("research/facets_qc.txt"))
        self.assertEqual(e[MetaDictMap.FIT_STATUS], "auto_qc_best_fit#default#auto-qc")
        tumor8, _ = self.info["twon"]["ultra_top"]
        self.assertEqual(m.master_file_dict[tumor8][MetaDictMap.FIT_STATUS], "not_reviewed#default#NA")  # ultra skipped
        self.assertEqual(m.build_counts["twon_samples"], 2)
        self.assertEqual(m.build_counts["skipped_ultra"], 3)

    def test_2n_fit_class_clinical_and_both(self):
        tumor, pair = self.info["twon"]["full"]
        m = self.meta(fit_class="clinical")
        self.assertEqual(m.master_file_dict[tumor][MetaDictMap.FIT_CLASS], "clinical")
        self.assertTrue(m.master_file_dict[tumor][MetaDictMap.SAMPLE_BASE_DIR].endswith("/clinical"))
        self.assertNotIn(tumor + "#clinical", m.master_file_dict)
        m = self.meta(fit_class="both")
        self.assertNotIn(tumor, m.master_file_dict)
        self.assertEqual(m.master_file_dict[tumor + "#clinical"][MetaDictMap.FIT_CLASS], "clinical")
        self.assertEqual(m.master_file_dict[tumor + "#research"][MetaDictMap.FIT_CLASS], "research")
        self.assertEqual(m.long_id_map[tumor], [pair])
        self.assertEqual(m.pairIdForKey(tumor + "#research"), pair)
        self.assertEqual(m.fit_map[tumor][0], "AutoQC_Best")       # short-id entry mirrors the first class
        # standard samples are untouched by fit_class
        std = self.info["standard"]["reviewed_best"][0]
        self.assertIn(std, m.master_file_dict)
        self.assertEqual(m.master_file_dict[std][MetaDictMap.FIT_CLASS], "standard")

    def test_all_fits_mode_keys_and_grouping(self):
        m = self.meta(single=False, fit_class="both")
        tumor, pair = self.info["twon"]["full"]
        keys = sorted(k for k in m.master_file_dict if k.startswith(pair))
        self.assertEqual(keys, [pair + "#clinical/alt_dipLogR_-0.1", pair + "#clinical/default",
                                pair + "#research/alt_dipLogR_-0.1", pair + "#research/default"])
        self.assertEqual(m.master_file_dict[pair + "#research/default"][MetaDictMap.FIT_CLASS], "research")
        self.assertEqual(m.master_file_dict[pair + "#research/default"][MetaDictMap.FIT_TIER], "autoqc")
        self.assertEqual(m.master_file_dict[pair + "#research/alt_dipLogR_-0.1"][MetaDictMap.FIT_TIER], "unreviewed")
        self.assertEqual(sorted(m.long_id_map[tumor]), keys)
        std_tumor, std_pair = self.info["standard"]["reviewed_best"]
        self.assertEqual(sorted(m.long_id_map[std_tumor]), [std_pair + "#alt_dipLogR_-0.1", std_pair + "#default"])
        self.assertEqual(m.fit_map[std_tumor][0], "Best")                                  # policy-preferred label
        self.assertEqual(m.fit_map[std_pair + "#default"][0], "Default")                   # per-run label
        self.assertEqual(m.pairIdForKey(std_pair + "#default"), std_pair)
        # all fits with ultra included: ultra dirs lack qc/gene files -> counted as missing, never crash
        m2 = self.meta(single=False, fit_class="both", policy=FitPolicy(include_ultra_fits=True), mode="hisens")
        self.assertEqual(m2.build_counts["skipped_ultra"], 0)
        self.assertEqual(m2.build_counts["missing_files"], 3)
        self.assertNotIn(pair + "#research/default/ultra", m2.master_file_dict)

    def test_duplicate_manifest_rows(self):
        tumor, pair = self.info["standard"]["duplicate_rows"]
        m = self.meta()
        self.assertEqual(m.master_file_dict[tumor][MetaDictMap.FIT_STATUS], "reviewed_best_fit#default#richara4/auto-qc-script")
        m = self.meta(single=False)
        keys = sorted(k for k in m.master_file_dict if k.startswith(pair))
        self.assertEqual(keys, [pair + "#alt_dipLogR_0.06", pair + "#default"])
        # the reviewed row wins over the later not_reviewed row for the same fit
        self.assertEqual(m.master_file_dict[pair + "#default"][MetaDictMap.FIT_STATUS], "reviewed_best_fit#default#richara4/auto-qc-script")
        self.assertEqual(m.master_file_dict[pair + "#default"][MetaDictMap.FIT_TIER], "autoqc")   # auto-qc-script reviewer
        self.assertEqual(m.fit_map[pair + "#default"][0], "AutoQC_Best")
        self.assertEqual(m.build_counts["duplicate_rows"], 1)
        self.assertEqual(m.build_counts["runs_registered"], len(m.master_file_dict))

    def test_auto_reviewer_classification_on_meta(self):
        tumor, pair = self.info["standard"]["duplicate_rows"]            # reviewed_best_fit by richara4/auto-qc-script
        human_tumor = self.info["standard"]["reviewed_best"][0]          # reviewed_best_fit by alice
        m = self.meta(policy=FitPolicy.reviewed_only(), defaults=False)
        self.assertIn(human_tumor, m.master_file_dict)
        self.assertNotIn(tumor, m.master_file_dict)
        m = self.meta(policy=FitPolicy.autoqc_only(), defaults=False)
        self.assertIn(tumor, m.master_file_dict)
        self.assertNotIn(human_tumor, m.master_file_dict)
        self.assertEqual(m.master_file_dict[tumor][MetaDictMap.FIT_STATUS], "reviewed_best_fit#default#richara4/auto-qc-script")
        self.assertEqual(m.fit_map[tumor][0], "AutoQC_Best")
        # a custom pattern can widen or narrow what counts as automated
        m = self.meta(policy=FitPolicy.reviewed_only(auto_reviewer_pattern=r"^$"), defaults=False)
        self.assertIn(tumor, m.master_file_dict)

    # ---- dataset -----------------------------------------------------------------------
    def test_dataset_single_run(self):
        m = self.meta(fit_class="both")
        ds = FacetsDataset(m)
        quiet(ds.buildFacetsDataset)
        self.assertEqual(len(ds.runList), 10)
        self.assertEqual(len(ds.sampleList), 8)
        tumor, pair = self.info["twon"]["full"]
        sample = ds.sampleList[pair]
        self.assertEqual(sorted(r.fit_class for r in sample.runs), ["clinical", "research"])
        run = [r for r in ds.runList if r.id == tumor + "#research"][0]
        self.assertEqual((run.fit_name, run.review_status, run.review_author, run.purity_hisens),
                         ("default", "auto_qc_best_fit", "auto-qc", "purity"))
        self.assertEqual(run.cancerType, m.cancer_type_map[tumor])
        self.assertTrue(run.facets_qc)
        self.assertEqual(run.cval, "100")
        self.assertEqual(float(run.ploidy), 2.1)
        self.assertGreater(len(run.segments), 0)
        self.assertGreater(len(run.genes), 0)
        calls = {g.gene: run.get_cn_call(g.tcn, g.lcn) for g in run.genes}
        self.assertEqual(calls["TP53"], "HETLOSS")
        self.assertEqual(calls["CDKN2A"], "HOMDEL")
        self.assertEqual(calls["EGFR"], "DIPLOID")
        std_tumor, std_pair = self.info["standard"]["reviewed_best"]
        self.assertIn(std_pair, ds.sampleList)
        self.assertEqual(ds.sampleList[std_pair].runs[0].id, std_tumor)

    def test_dataset_filters_in_all_fits_mode(self):
        m = self.meta(single=False, fit_class="both")
        ds = FacetsDataset(m)
        ds.setCancerTypeDetailFilter(["Lung Adenocarcinoma"])      # used to crash on an attribute typo
        quiet(ds.buildFacetsDataset)
        self.assertGreater(len(ds.runList), 0)
        self.assertTrue(all(r.cancerTypeDetail == "Lung Adenocarcinoma" for r in ds.runList))
        ds = FacetsDataset(m)
        ds.setClinicalPurityFilter(25, 100)                           # NA purities are excluded, not crashed on
        quiet(ds.buildFacetsDataset)
        self.assertTrue(all(float(r.clinicalPurity) >= 25 for r in ds.runList))
        ds = FacetsDataset(m)
        ds.setPurityFilter(0.5, 0.7)
        ds.setFacetsQCFilter(True)
        quiet(ds.buildFacetsDataset)
        self.assertGreater(len(ds.runList), 0)
        no_clin = self.meta(clinical=False)
        with self.assertRaises(FacetsConfigError):
            FacetsDataset(no_clin).setCancerTypeFilter(["Breast Cancer"])

    def test_hisens_mode_reads_hisens_cval(self):
        m = self.meta(mode="hisens")
        ds = FacetsDataset(m)
        quiet(ds.buildFacetsDataset)
        self.assertTrue(all(r.cval == "50" and r.purity_hisens == "hisens" for r in ds.runList))
        self.assertGreater(len(ds.runList), 0)

    def test_writers(self):
        m = self.meta(fit_class="both")
        ds = FacetsDataset(m)
        quiet(ds.buildFacetsDataset)
        out = os.path.join(self.tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        quiet(ds.writeReport, os.path.join(out, "report.tsv"))
        lines = open(os.path.join(out, "report.tsv")).read().splitlines()
        self.assertEqual(len(lines), 1 + len(ds.runList))
        self.assertIn("AutoQC_Best", "\n".join(lines))
        ds.writePurityCFs(os.path.join(out, "cf.tsv"))
        self.assertEqual(len(open(os.path.join(out, "cf.tsv")).read().splitlines()), 1 + len(ds.runList))
        ds.writeDatasetSummary(outdir=out)
        self.assertTrue(any(f.startswith("FacetsDataset_summary_") for f in os.listdir(out)))
        p = ds.createHistogram("purity", outdir=out)
        self.assertTrue(os.path.isfile(p))
        p = ds.createHistogram("fit_class", outdir=out)
        self.assertTrue(os.path.isfile(p))
        with self.assertRaises(FacetsConfigError):
            ds.createHistogram("bogus", outdir=out)

    def test_sample_list_file_filters_by_tumor_id(self):
        tumor, pair = self.info["standard"]["reviewed_best"]
        lst = os.path.join(self.tmp.name, "samples.lst")
        with open(lst, "w") as fh:
            fh.write("# a comment\n%s\n\n%s\n" % (tumor, self.info["twon"]["full"][1]))   # tumor id and a pair id
        m = FacetsMeta(self.clin, self.repo, "purity")
        m.setSingleRunPerSample(True, True)
        m.selectSamplesFromFile(lst)
        quiet(m.buildFacetsMeta)
        self.assertEqual(sorted(m.master_file_dict), sorted([tumor, self.info["twon"]["full"][0]]))
        with self.assertRaises(FacetsRepositoryError):
            m.selectSamplesFromFile(os.path.join(self.tmp.name, "nope.lst"))

    def test_persistence_roundtrip_and_old_format_warning(self):
        store = os.path.join(self.tmp.name, "meta.dat")
        m = FacetsMeta(self.clin, self.repo, "purity", persist_data=store)
        m.setSingleRunPerSample(True, True)
        quiet(m.buildFacetsMeta)
        m2 = FacetsMeta(self.clin, self.repo, "purity", persist_data=store)
        m2.setSingleRunPerSample(True, True)
        quiet(m2.buildFacetsMeta)
        self.assertEqual(m2.master_file_dict, m.master_file_dict)
        self.assertTrue(m2.selectSingleRun)
        # simulate an old (0.5) pickle: no _persist_version
        old = os.path.join(self.tmp.name, "old.dat")
        d = dict(m.__dict__)
        d.pop("_persist_version")
        o = _OldPickledMeta(); o.__dict__.update(d)
        with open(old, "wb") as fh:
            pickle.dump(o, fh)
        m3 = FacetsMeta(self.clin, self.repo, "hisens", persist_data=old)
        m3.setSingleRunPerSample(True, False)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            m3.buildFacetsMeta()
        self.assertIn("older facetsAPI", buf.getvalue())
        self.assertEqual(m3.hisens_vs_purity, "hisens")                 # constructor value kept
        self.assertEqual(m3.fit_policy.tiers, ("reviewed", "autoqc"))   # setter value kept
        self.assertEqual(m3.master_file_dict, m.master_file_dict)


if __name__ == "__main__":
    unittest.main()

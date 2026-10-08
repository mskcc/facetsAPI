import os
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

from facetsAPI.exceptions import FacetsConfigError  # noqa: E402
from facetsAPI.layout import read_manifest  # noqa: E402
from facetsAPI.fitpolicy import FitPolicy, tier_of, label_of  # noqa: E402


class FitPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.repo = cls.info["repo"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def manifest(self, group, name, cls=None):
        tumor, pair = self.info[group][name][0:2]
        d = os.path.join(self.repo, tumor[0:7], pair)
        if cls:
            d = os.path.join(d, cls)
        return read_manifest(os.path.join(d, "facets_review.manifest")), d

    def names(self, cands):
        return [(c.fit_name, c.tier, c.label) for c in cands]

    def test_default_policy_order(self):
        df, _ = self.manifest("standard", "reviewed_best")
        self.assertEqual(self.names(FitPolicy().candidates(df)),
                         [("alt_dipLogR_-0.1", "reviewed", "Best"), ("default", "default", "Default")])
        df, _ = self.manifest("standard", "autoqc_only")
        self.assertEqual(self.names(FitPolicy().candidates(df)),
                         [("alt_dipLogR_0.2", "autoqc", "AutoQC_Best"), ("default", "autoqc", "AutoQC_Pass")])
        df, _ = self.manifest("standard", "default_only")
        self.assertEqual(self.names(FitPolicy().candidates(df)), [("default", "default", "Default")])

    def test_only_policies(self):
        a, _ = self.manifest("standard", "reviewed_best")
        b, _ = self.manifest("standard", "autoqc_only")
        c, _ = self.manifest("standard", "default_only")
        self.assertEqual(self.names(FitPolicy.reviewed_only().candidates(b)), [])
        self.assertEqual(self.names(FitPolicy.reviewed_only().candidates(c)), [])
        self.assertEqual(self.names(FitPolicy.autoqc_only().candidates(a)), [])
        self.assertEqual(self.names(FitPolicy.autoqc_only().candidates(b)),
                         [("alt_dipLogR_0.2", "autoqc", "AutoQC_Best"), ("default", "autoqc", "AutoQC_Pass")])
        self.assertEqual(self.names(FitPolicy.default_only().candidates(a)), [("default", "default", "Default")])
        self.assertEqual(self.names(FitPolicy.default_only().candidates(b)), [("default", "default", "AutoQC_Pass")])  # tier = where the policy found it; label = its status
        self.assertEqual(self.names(FitPolicy.reviewed_then_autoqc().candidates(c)), [])
        self.assertFalse(FitPolicy.reviewed_then_autoqc().allows_defaults)

    def test_ultra_excluded_by_default_and_fallthrough(self):
        df, research = self.manifest("twon", "ultra_top", "research")
        self.assertEqual(self.names(FitPolicy().candidates(df)), [("default", "default", "Default")])
        self.assertEqual(self.names(FitPolicy.autoqc_only().candidates(df)), [])
        inc = FitPolicy(include_ultra_fits=True)
        self.assertEqual(self.names(inc.candidates(df)),
                         [("alt_dipLogR_0.3/ultra", "autoqc", "AutoQC_Pass"), ("default", "default", "Default")])
        cand, files, tried = inc.select_fit(df, "purity", research)
        self.assertEqual(cand.fit_name, "default")          # ultra has no purity files -> fell through
        self.assertEqual([t.fit_name for t in tried], ["alt_dipLogR_0.3/ultra"])
        self.assertTrue(files["out"].endswith("_purity.out"))
        cand, files, tried = FitPolicy.autoqc_only(include_ultra_fits=True).select_fit(df, "hisens", research)
        self.assertIsNone(cand)                              # ultra dir lacks qc/gene files
        self.assertEqual(len(tried), 1)

    def test_2n_full_research_and_not_selected(self):
        df, research = self.manifest("twon", "full", "research")
        cands = FitPolicy(include_ultra_fits=True).candidates(df)
        self.assertEqual([c.fit_name for c in cands], ["default", "alt_dipLogR_-0.1/ultra"])
        self.assertNotIn("Not selected", [c.fit_name for c in cands])
        cand, files, tried = FitPolicy().select_fit(df, "purity", research)
        self.assertEqual((cand.fit_name, cand.label, cand.reviewed_by), ("default", "AutoQC_Best", "auto-qc"))
        self.assertEqual(cand.fit_status, "auto_qc_best_fit#default#auto-qc")
        self.assertEqual(tried, [])

    def test_no_fit_cases(self):
        df, d = self.manifest("standard", "default_only")
        cand, files, tried = FitPolicy.reviewed_only().select_fit(df, "purity", d)
        self.assertIsNone(cand)
        self.assertEqual(tried, [])
        df, d = self.manifest("twon", "ultra_top", "clinical")   # default has auto_qc_no_fit status
        cand, _, _ = FitPolicy().select_fit(df, "purity", d)
        self.assertEqual((cand.fit_name, cand.tier), ("default", "default"))

    def test_with_defaults_and_validation(self):
        p = FitPolicy.reviewed_only()
        self.assertEqual(p.with_defaults(True).tiers, ("reviewed", "default"))
        self.assertEqual(FitPolicy().with_defaults(False).tiers, ("reviewed", "autoqc"))
        self.assertEqual(FitPolicy.default_only().with_defaults(False).tiers, ("default",))
        with self.assertRaises(FacetsConfigError):
            FitPolicy(tiers=("reviewed", "bogus"))
        with self.assertRaises(FacetsConfigError):
            FitPolicy(tiers=())
        with self.assertRaises(FacetsConfigError):
            FitPolicy(tie_break="random")
        self.assertEqual(FitPolicy(tiers="autoqc").tiers, ("autoqc",))

    def test_tier_and_label(self):
        self.assertEqual(tier_of("reviewed_best_fit", "alt"), "reviewed")
        self.assertEqual(tier_of("auto_qc_pass", "alt"), "autoqc")
        self.assertEqual(tier_of("not_reviewed", "default"), "default")
        self.assertEqual(tier_of("auto_qc_no_fit", "alt"), "unreviewed")
        self.assertEqual(label_of("reviewed_acceptable_fit", "x"), "Acceptable")
        self.assertEqual(label_of("not_reviewed", "default"), "Default")
        self.assertEqual(label_of("not_reviewed", "alt"), "Unreviewed")

    def test_representative_rows(self):
        df = pd.DataFrame({
            "path": ["/p/"] * 5, "fit_name": ["default", "alt", "default", "refit", "refit"],
            "review_status": ["reviewed_best_fit", "not_reviewed", "not_reviewed", "reviewed_best_fit", "reviewed_best_fit"],
            "reviewed_by": ["auto-qc-script", "NA", "NA", "alice", "bob"],
            "date_reviewed": ["2022-12-15 02:03:19", "NA", "NA", "2020-06-01 16:09:19", "2024-01-01 00:00:00"],
        })
        rows, dropped = FitPolicy().representative_rows(df)
        self.assertEqual(dropped, 2)
        self.assertEqual(list(rows["fit_name"]), ["default", "alt", "refit"])
        self.assertEqual(rows.loc[rows["fit_name"] == "default", "review_status"].iloc[0], "reviewed_best_fit")
        self.assertEqual(rows.loc[rows["fit_name"] == "refit", "reviewed_by"].iloc[0], "bob")   # newest review wins

    def test_newest_reviewed_tie_break(self):
        df = pd.DataFrame({
            "path": ["/p/"] * 3, "fit_name": ["old", "new", "default"],
            "review_status": ["reviewed_best_fit", "reviewed_best_fit", "not_reviewed"],
            "reviewed_by": ["a", "b", "NA"], "date_reviewed": ["2024-01-01", "2025-06-01", "NA"],
        })
        self.assertEqual([c.fit_name for c in FitPolicy().candidates(df)], ["old", "new", "default"])
        self.assertEqual([c.fit_name for c in FitPolicy(tie_break="newest_reviewed").candidates(df)],
                         ["new", "old", "default"])


if __name__ == "__main__":
    unittest.main()

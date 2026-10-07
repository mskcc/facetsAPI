import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

from facetsAPI.exceptions import FacetsConfigError, ManifestError  # noqa: E402
from facetsAPI import layout  # noqa: E402


class LayoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.repo = cls.info["repo"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def sdir(self, tumor, pair):
        return os.path.join(self.repo, tumor[0:7], pair)

    def test_discover_standard(self):
        tumor, pair = self.info["standard"]["reviewed_best"]
        classes = layout.discover_sample(self.sdir(tumor, pair))
        self.assertEqual([c.name for c in classes], ["standard"])
        self.assertFalse(classes[0].is_2n)
        self.assertTrue(os.path.isfile(classes[0].manifest_path))
        self.assertTrue(os.path.isfile(classes[0].facets_qc_path))

    def test_discover_2n(self):
        tumor, pair = self.info["twon"]["full"]
        classes = layout.discover_sample(self.sdir(tumor, pair))
        self.assertEqual([c.name for c in classes], ["clinical", "research"])
        self.assertTrue(all(c.is_2n for c in classes))
        self.assertTrue(classes[1].class_dir.endswith("/research"))

    def test_discover_none(self):
        tumor, pair = self.info["standard"]["no_manifest"]
        self.assertEqual(layout.discover_sample(self.sdir(tumor, pair)), [])

    def test_select_classes(self):
        tumor, pair = self.info["twon"]["full"]
        classes = layout.discover_sample(self.sdir(tumor, pair))
        self.assertEqual([c.name for c in layout.select_classes(classes, "research")], ["research"])
        self.assertEqual([c.name for c in layout.select_classes(classes, "clinical")], ["clinical"])
        self.assertEqual([c.name for c in layout.select_classes(classes, "both")], ["clinical", "research"])
        tumor, pair = self.info["standard"]["reviewed_best"]
        std = layout.discover_sample(self.sdir(tumor, pair))
        self.assertEqual(layout.select_classes(std, "clinical"), std)  # standard unaffected

    def test_validate_fit_class(self):
        with self.assertRaises(FacetsConfigError):
            layout.validate_fit_class("ultra")
        self.assertEqual(layout.validate_fit_class("both"), "both")
        with self.assertRaises(FacetsConfigError):
            layout.validate_mode("hisense")

    def test_read_manifest_with_and_without_comment(self):
        tumor, pair = self.info["standard"]["reviewed_best"]
        df = layout.read_manifest(os.path.join(self.sdir(tumor, pair), "facets_review.manifest"))
        self.assertEqual(sorted(df["fit_name"]), ["alt_dipLogR_-0.1", "default"])
        tumor, pair = self.info["standard"]["no_comment_acceptable"]
        df2 = layout.read_manifest(os.path.join(self.sdir(tumor, pair), "facets_review.manifest"))
        self.assertEqual(list(df2["fit_name"]), ["default"])
        self.assertEqual(list(df2["review_status"]), ["reviewed_acceptable_fit"])
        self.assertEqual(list(df2["reviewed_by"]), ["carol"])

    def test_read_manifest_na_reviewer_and_missing_column(self):
        tumor, pair = self.info["standard"]["default_only"]
        df = layout.read_manifest(os.path.join(self.sdir(tumor, pair), "facets_review.manifest"))
        self.assertEqual(df.loc[df["fit_name"] == "default", "reviewed_by"].iloc[0], "NA")
        bad = os.path.join(self.tmp.name, "bad.manifest")
        with open(bad, "w") as fh:
            fh.write("# comment\n# another\nsample\tpath\treview_status\n")
            fh.write("x\t/p/\tnot_reviewed\n")
        with self.assertRaises(ManifestError) as cm:
            layout.read_manifest(bad)
        self.assertIn("fit_name", str(cm.exception))
        with self.assertRaises(ManifestError):
            layout.read_manifest(os.path.join(self.tmp.name, "does_not_exist.manifest"))

    def test_ids(self):
        self.assertEqual(layout.short_sample_id("P-0000001-T01-IM6"), "P-0000001-T01-IM6")
        self.assertEqual(layout.short_sample_id("P-0000001-T01-IM6#research"), "P-0000001-T01-IM6")
        self.assertEqual(layout.short_sample_id("P-0000001-T01-IM6_P-0000001-N01-IM6"), "P-0000001-T01-IM6")
        self.assertEqual(layout.short_sample_id("P-0000001-T01-IM6_P-0000001-N01-IM6#research/default/ultra"),
                         "P-0000001-T01-IM6")
        self.assertEqual(layout.pair_id("P-0000001-T01-IM6_P-0000001-N01-IM6#research/default"),
                         "P-0000001-T01-IM6_P-0000001-N01-IM6")
        self.assertTrue(layout.is_ultra("default/ultra"))
        self.assertTrue(layout.is_ultra("alt_dipLogR_-0.1/ultra"))
        self.assertFalse(layout.is_ultra("default"))
        self.assertEqual(layout.fit_dir("/x/research", "default/ultra"), "/x/research/default/ultra")

    def test_locate_fit_files(self):
        tumor, pair = self.info["twon"]["full"]
        research = os.path.join(self.sdir(tumor, pair), "research")
        files, missing = layout.locate_fit_files(layout.fit_dir(research, "default"), "purity")
        self.assertEqual(missing, [])
        self.assertTrue(files["out"].endswith(pair + "_purity.out"))
        self.assertTrue(files["qc"].endswith(".qc.txt"))
        self.assertTrue(files["nonsignedout"].endswith("nonsignedout.ccf.maf"))
        self.assertTrue(files["ccf"].endswith(pair + ".ccf.maf"))
        files, missing = layout.locate_fit_files(layout.fit_dir(research, "default/ultra"), "purity")
        self.assertIsNone(files)
        self.assertIn("out", missing)
        files, missing = layout.locate_fit_files(layout.fit_dir(research, "default/ultra"), "hisens")
        self.assertIsNone(files)                      # ultra dirs carry no qc/gene files on disk
        self.assertEqual(sorted(missing), ["gene", "qc"])
        tumor, pair = self.info["standard"]["no_comment_acceptable"]
        files, missing = layout.locate_fit_files(layout.fit_dir(self.sdir(tumor, pair), "default"), "purity")
        self.assertEqual(missing, [])
        self.assertEqual(files["rdata"], "")          # Rdata is optional


if __name__ == "__main__":
    unittest.main()

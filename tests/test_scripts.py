"""Run the salvaged analysis scripts end to end on the synthetic repository."""
import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import unittest
import unittest.mock

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class ScriptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.samples = [cls.info["standard"]["reviewed_best"][0], cls.info["standard"]["autoqc_only"][0],
                       cls.info["twon"]["full"][0]]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_get_state_percents(self):
        gsp = load("get_state_percents")
        self.assertEqual(gsp.classify_hetloss_homdel(1, 0), "HETLOSS")
        self.assertEqual(gsp.classify_hetloss_homdel(1, -1), "HETLOSS")
        self.assertEqual(gsp.classify_hetloss_homdel(0, 0), "HOMDEL")
        self.assertIsNone(gsp.classify_hetloss_homdel(2, 1))
        out = os.path.join(self.tmp.name, "states.tsv")
        lst = os.path.join(self.tmp.name, "samples.lst")
        with open(lst, "w") as fh:
            fh.write("\n".join(self.samples) + "\n")
        argv = ["get_state_percents", "--facets-dir", self.info["repo"], "--samples-file", lst, "-o", out]
        with unittest.mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
            gsp.main()
        rows = [ln.split("\t") for ln in open(out).read().splitlines()]
        self.assertEqual(rows[0], ["sample", "pct_HETLOSS", "pct_HOMDEL", "bp_HETLOSS", "bp_HOMDEL", "total_bp"])
        self.assertEqual(len(rows), 1 + len(self.samples))
        for r in rows[1:]:
            self.assertGreater(float(r[1]), 0.0)   # fixture has (1,0) segments on chr1 and chr17
            self.assertGreater(float(r[2]), 0.0)   # and a (0,0) segment on chr9
            self.assertGreater(int(r[5]), 0)

    def test_cin_scores_cli(self):
        cin = load("cin_scores")
        self.assertEqual(cin.arm_sort_key("2q"), (2, 1))
        self.assertEqual(cin.arm_sort_key("Xp"), (23, 0))
        prefix = os.path.join(self.tmp.name, "cin")
        argv = ["--facets-dir", self.info["repo"], "--clinical-sample-file", self.info["clinical_file"],
                "--sample", self.samples[0], "--sample", self.samples[1], "-o", prefix]
        with contextlib.redirect_stdout(io.StringIO()):
            summary = cin.main(argv)
        for suffix in ("_cin_full.tsv", "_cin_arm.tsv", "_candidate_arms.tsv"):
            self.assertTrue(os.path.isfile(prefix + suffix), suffix)
        full = open(prefix + "_cin_full.tsv").read().splitlines()
        self.assertEqual(len(full), 1 + 2)
        self.assertIn("score", full[0])
        arm = open(prefix + "_cin_arm.tsv").read().splitlines()
        self.assertGreater(len(arm), 3)
        self.assertEqual(sorted(summary), sorted(self.samples[0:2]))
        cands = open(prefix + "_candidate_arms.tsv").read().splitlines()
        self.assertEqual(cands[0].split("\t")[0:4], ["sample", "arm", "score", "breakpoints"])


if __name__ == "__main__":
    unittest.main()

import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

from facetsAPI import FacetsMeta, FacetsDataset, ExtTools, FPTools, FacetsConfig, MetaDictMap  # noqa: E402
from facetsAPI.exceptions import ExternalToolError, FacetsConfigError  # noqa: E402

IRIS_PREFIX = "module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3"


def completed(rc=0, out="", err=""):
    return subprocess.CompletedProcess(args=["bash", "-lc", "x"], returncode=rc, stdout=out, stderr=err)


class FacetsConfigTests(unittest.TestCase):
    def test_defaults_point_into_repo(self):
        c = FacetsConfig()
        self.assertEqual(c.rscript, "Rscript")
        self.assertEqual(c.r_prefix, "")
        self.assertTrue(c.fp_wrapper_script.endswith("R/fp_generate_genomic_anno_wrapper.R"))
        self.assertTrue(os.path.isfile(c.fp_wrapper_script))
        self.assertTrue(os.path.isfile(os.path.join(c.ascets_dir, "run_ascets.R")))
        self.assertIsNone(c.fp_config_file)

    def test_env_precedence(self):
        env = {"FACETSAPI_R_PREFIX": IRIS_PREFIX, "FACETSAPI_RSCRIPT": "/opt/R/bin/Rscript", "FACETSAPI_FP_CONFIG": "/cfg.json"}
        c = FacetsConfig.from_env(environ=env)
        self.assertEqual((c.r_prefix, c.rscript, c.fp_config_file), (IRIS_PREFIX, "/opt/R/bin/Rscript", "/cfg.json"))
        c = FacetsConfig.from_env(environ=env, rscript="Rscript", fp_config_file=None)
        self.assertEqual(c.rscript, "Rscript")              # explicit wins
        self.assertEqual(c.fp_config_file, "/cfg.json")     # None override keeps env
        with mock.patch.dict(os.environ, {"FACETSAPI_R_LIB_LOC": "/libs"}, clear=False):
            self.assertEqual(FacetsConfig.from_env().r_lib_loc, "/libs")
        with self.assertRaises(FacetsConfigError):
            FacetsConfig.from_env(bogus="x")

    def test_compose_command(self):
        c = FacetsConfig()
        self.assertEqual(c.compose_r_command(["/a/b.R", "-i", "x y"]), "Rscript /a/b.R -i 'x y'")
        c = FacetsConfig(r_prefix=IRIS_PREFIX, rscript="Rscript")
        self.assertEqual(c.compose_r_command(["/a/b.R"]), IRIS_PREFIX + " && Rscript /a/b.R")
        with self.assertRaises(FacetsConfigError):
            FacetsConfig(rscript="").compose_r_command(["x"])

    def test_run_r_uses_bash_lc_and_raises_on_failure(self):
        c = FacetsConfig(r_prefix=IRIS_PREFIX)
        with mock.patch("facetsAPI.config.subprocess.run", return_value=completed(0, "ok")) as run:
            res = c.run_r(["/w.R", "a"])
            self.assertEqual(res.stdout, "ok")
            argv = run.call_args[0][0]
            self.assertEqual(argv[0:2], ["bash", "-lc"])
            self.assertEqual(argv[2], IRIS_PREFIX + " && Rscript /w.R a")
            self.assertTrue(run.call_args[1]["capture_output"])
        with mock.patch("facetsAPI.config.subprocess.run", return_value=completed(2, "", "Error in library(facetsPreview)")):
            with self.assertRaises(ExternalToolError) as cm:
                c.run_r(["/w.R"])
            self.assertEqual(cm.exception.returncode, 2)
            self.assertIn("library(facetsPreview)", str(cm.exception))
        with mock.patch("facetsAPI.config.subprocess.run", side_effect=FileNotFoundError("bash")):
            with self.assertRaises(ExternalToolError):
                c.run_r(["/w.R"])


class ToolsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.meta = FacetsMeta(cls.info["clinical_file"], cls.info["repo"], "purity", fit_class="both")
        cls.meta.setSingleRunPerSample(True, True)
        with contextlib.redirect_stdout(io.StringIO()):
            cls.meta.buildFacetsMeta()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_exttools_validation(self):
        with self.assertRaises(FacetsConfigError):
            ExtTools("not a meta")
        with self.assertRaises(FacetsConfigError):
            ExtTools(self.meta, config="nope")
        with self.assertRaises(NotImplementedError):
            ExtTools(self.meta).runFacets()
        with self.assertRaises(FacetsConfigError):
            ExtTools(self.meta).runAscets("")
        with self.assertRaises(FacetsConfigError):
            ExtTools(self.meta).runAscets(self.tmp.name, ref_genome_coords="hg18")

    def test_make_merged_file(self):
        out = os.path.join(self.tmp.name, "merged.seg")
        ExtTools(self.meta).makeMergedFile(MetaDictMap.ADJUSTED_SEG_FILE, out)
        lines = open(out).read().splitlines()
        self.assertEqual(lines[0], "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean")
        self.assertEqual(sum(1 for l in lines if l.startswith("ID\t")), 1)   # one header
        n_runs = len(self.meta.master_file_dict)
        # 2n 'both' loads clinical+research for the same pair: identical basenames are de-duplicated
        n_unique = len({os.path.basename(v[MetaDictMap.ADJUSTED_SEG_FILE]) for v in self.meta.master_file_dict.values()})
        self.assertLess(n_unique, n_runs)
        self.assertEqual(len(lines) - 1, n_unique * 7)

    def test_run_ascets_composes_and_logs(self):
        out_dir = os.path.join(self.tmp.name, "ascets_out")
        cfg = FacetsConfig(r_prefix=IRIS_PREFIX)
        tools = ExtTools(self.meta, config=cfg)
        with mock.patch("facetsAPI.config.subprocess.run", return_value=completed(0, "ran", "warn")) as run, \
                contextlib.redirect_stdout(io.StringIO()):
            prefix = tools.runAscets(out_dir, ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False,
                                     arm_alt_frac_thresh=0.7, use_adjusted=False)
        cmd = run.call_args[0][0][2]
        self.assertTrue(cmd.startswith(IRIS_PREFIX + " && Rscript "))
        self.assertIn("run_ascets.R -i " + os.path.join(out_dir, "merged_seg_data.seg"), cmd)
        self.assertIn("genomic_arm_coordinates_hg19.txt -m 0.5 -k F -a 0.7 -o FacetsAscetsRun_hg19_breadth0.5_filterNoise_armAltFracThresh0.7_facetsUnadjustedSeg_", cmd)
        self.assertTrue(prefix.startswith("FacetsAscetsRun_hg19"))
        self.assertEqual(open(os.path.join(out_dir, "runLog.out")).read(), "ran")
        self.assertEqual(open(os.path.join(out_dir, "runLog.err")).read(), "warn")
        self.assertTrue(os.path.isfile(os.path.join(out_dir, "merged_seg_data.seg")))
        with mock.patch("facetsAPI.config.subprocess.run", return_value=completed(1, "", "boom")), \
                contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ExternalToolError):
                tools.runAscets(out_dir)

    def test_fptools(self):
        with self.assertRaises(FacetsConfigError):
            FPTools().runGenerateGenomicAnnotations("T", "T_N", "/p")       # no config file anywhere
        fp = FPTools("/cfg/facets_preview.json", "/m/cbio.maf", "/m/cbio.nonsigned.maf",
                     config=FacetsConfig(r_prefix=IRIS_PREFIX, r_lib_loc="/R_libs"))
        self.assertEqual(fp.config_file, "/cfg/facets_preview.json")       # legacy attribute
        with mock.patch("facetsAPI.config.subprocess.run", return_value=completed(0)) as run, \
                contextlib.redirect_stdout(io.StringIO()):
            fp.runGenerateGenomicAnnotations("P-0000001-T01-IM6", "P-0000001-T01-IM6_P-0000001-N01-IM6", "/s/path")
        cmd = run.call_args[0][0][2]
        self.assertTrue(cmd.startswith(IRIS_PREFIX + " && Rscript "))
        self.assertTrue(cmd.endswith("fp_generate_genomic_anno_wrapper.R /cfg/facets_preview.json P-0000001-T01-IM6 "
                                     "P-0000001-T01-IM6_P-0000001-N01-IM6 /s/path /m/cbio.maf /m/cbio.nonsigned.maf /R_libs"), cmd)
        with mock.patch.dict(os.environ, {"FACETSAPI_FP_CONFIG": "/env/cfg.json"}):
            self.assertEqual(FPTools().config_file, "/env/cfg.json")
        bad = FPTools("/cfg.json", config=FacetsConfig(fp_wrapper_script="/nope.R"))
        with self.assertRaises(FacetsConfigError):
            bad.runGenerateGenomicAnnotations("a", "b", "c")

    def test_backups(self):
        with contextlib.redirect_stdout(io.StringIO()):
            n = FPTools("/cfg.json").backupManifests(self.meta)
            q = FPTools("/cfg.json").backupFacetsQC(self.meta)
            ds = FacetsDataset(self.meta)
            ds.buildFacetsDataset()
            m = ds.backupMafs(self.meta)
        manifests = {v[MetaDictMap.MANIFEST_FILE] for v in self.meta.master_file_dict.values()}
        self.assertEqual(n, len(manifests))
        self.assertEqual(q, len(manifests))
        self.assertGreater(m, 0)
        for man in manifests:
            self.assertTrue(any(f.startswith("facets_review.manifest.") and f.endswith(".bak")
                                for f in os.listdir(os.path.dirname(man))))
        with self.assertRaises(FacetsConfigError):
            FPTools("/cfg.json").backupManifests("x")


if __name__ == "__main__":
    unittest.main()

"""Request spec validation and script generation, including rendered scripts run on the synthetic repository."""
import copy
import csv
import glob
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from make_fixture_repo import build  # noqa: E402

from facetsAPI import spec as S  # noqa: E402

CHECKOUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sid(number):
    """A synthetic tumor sample ID, built at runtime so no real-looking IDs are committed."""
    return "-".join(["P", str(number).zfill(7), "T01", "IM6"])


def r_has_ascets_packages():
    code = 'q(status = !all(sapply(c("dplyr", "tidyr", "tibble", "data.table"), requireNamespace, quietly = TRUE)))'
    try:
        return subprocess.run(["Rscript", "-e", code], capture_output=True).returncode == 0
    except OSError:
        return False


def base_spec(**overrides):
    spec = {
        "spec_version": S.SPEC_VERSION,
        "name": "test_job",
        "source": {"facets_dir": "/repo/all/", "clinical_file": "/repo/clinical.txt", "mode": "purity",
                   "fit_class": "research"},
        "fits": {"single_run": True, "tiers": ["reviewed", "autoqc", "default"]},
        "samples": {"source": "all"},
        "filters": {},
        "outputs": [{"type": "report"}],
    }
    spec.update(overrides)
    return spec


def ctx(**overrides):
    values = {"job_id": "job-1", "user": "someone", "facetsapi_home": CHECKOUT, "facetsapi_sha": "abc123",
              "output_dir": "/out/test_job", "job_dir": "/jobs/job-1"}
    values.update(overrides)
    return values


def every_option_spec():
    return base_spec(
        samples={"source": "list", "ids": [sid(1)]},
        filters={
            "purity": {"min": 0.3}, "clinical_purity": {"min": 10, "max": 90}, "ploidy": {"min": 1, "max": 5},
            "diplogr": {"max": 2}, "cval": {"min": 1}, "fga": {"max": 0.9}, "frac_loh": {"min": 0},
            "tmb": {"min": 0}, "msi": {"max": 100}, "wgd": False, "facets_qc": True,
            "cancer_type": ["Breast Cancer"], "cancer_type_detailed": ["Invasive Ductal Carcinoma"],
            "oncotree": {"presets": ["breast_carcinoma"], "codes": ["LUAD"]},
        },
        outputs=[
            {"type": "report"},
            {"type": "gene_level", "genes": ["TP53"]},
            {"type": "alteration", "thresholds": {"alteration_arm_num": 10, "percent_arm_loh": 40.0}},
            {"type": "summary"},
            {"type": "histograms", "variables": ["purity", "fit_name"], "bins": 5},
            {"type": "purity_cf"},
            {"type": "merged", "files": ["adjusted_seg", "cncf"]},
            {"type": "cin", "levels": ["full", "arm", "chrom"], "gene_or_seg": "seg", "split_arms": True},
            {"type": "chromothripsis", "gene_or_seg": "seg", "min_breakpoints": 3, "min_ris_delta_pct": 0.5,
             "breakpoint_offset": 1, "bg_metric": "mean"},
            {"type": "state_percents", "include_x": True, "include_y": False},
            {"type": "ascets", "genome": "hg19", "min_arm_breadth": 0.5, "arm_alt_frac_thresh": 0.7,
             "keep_noise": False, "use_adjusted": True},
            {"type": "copy_folders"},
        ],
    )


class ValidateTests(unittest.TestCase):
    def assertInvalid(self, spec, fragment):
        errors, _ = S.validate(spec)
        self.assertTrue(any(fragment in e for e in errors), "%r not in %r" % (fragment, errors))

    def test_valid_specs(self):
        self.assertEqual(S.validate(base_spec())[0], [])
        self.assertEqual(S.validate(every_option_spec())[0], [])

    def test_not_an_object(self):
        self.assertInvalid([], "must be an object")

    def test_version_and_name(self):
        self.assertInvalid(base_spec(spec_version=99), "spec_version")
        self.assertInvalid(base_spec(name="../escape"), "Job name")

    def test_enums(self):
        spec = base_spec()
        spec["source"]["mode"] = "other"
        self.assertInvalid(spec, "source.mode")
        spec = base_spec()
        spec["source"]["fit_class"] = "other"
        self.assertInvalid(spec, "source.fit_class")
        self.assertInvalid(base_spec(fits={"single_run": True, "tiers": ["bogus"]}), "fits.tiers")
        self.assertInvalid(base_spec(fits={"single_run": True, "tiers": []}), "fits.tiers")
        self.assertInvalid(base_spec(fits={"single_run": True, "tiers": ["default", "default"]}), "fits.tiers")
        self.assertInvalid(base_spec(outputs=[{"type": "backupMafs"}]), "Unknown output")
        self.assertInvalid(base_spec(filters={"bogus": 1}), "Unknown filter")

    def test_unhashable_values_are_errors_not_crashes(self):
        self.assertInvalid(base_spec(outputs=[{"type": "merged", "files": [["cncf"]]}]), "Merged files")
        self.assertInvalid(base_spec(outputs=[{"type": "histograms", "variables": [{}]}]), "Histograms")
        self.assertInvalid(base_spec(outputs=[{"type": "cin", "levels": [[]]}]), "CIN")
        self.assertInvalid(base_spec(outputs=[{"type": "cin", "levels": ["full"], "gene_or_seg": []}]), "CIN")
        self.assertInvalid(base_spec(fits={"single_run": True, "tiers": ["default"], "tie_break": []}), "tie_break")
        self.assertInvalid(base_spec(outputs=[{"type": []}]), "Unknown output")

    def test_new_output_options(self):
        self.assertInvalid(base_spec(outputs=[{"type": "cin", "levels": []}]), "CIN")
        self.assertInvalid(base_spec(outputs=[{"type": "chromothripsis", "min_breakpoints": 1.5}]), "min_breakpoints")
        self.assertInvalid(base_spec(outputs=[{"type": "chromothripsis", "bg_metric": "mode"}]), "background")
        self.assertInvalid(base_spec(outputs=[{"type": "ascets", "genome": "mm10"}]), "genome")
        self.assertInvalid(base_spec(outputs=[{"type": "ascets", "min_arm_breadth": 2}]), "min_arm_breadth")
        self.assertInvalid(base_spec(outputs=[{"type": "state_percents", "include_x": "yes"}]), "include_x")

    def test_sample_scores_need_one_fit_per_sample(self):
        for output in ({"type": "cin", "levels": ["full"]}, {"type": "chromothripsis"}, {"type": "state_percents"}, {"type": "ascets"}):
            self.assertInvalid(base_spec(fits={"single_run": False, "tiers": ["default"]}, outputs=[output]), "one fit per sample")
            spec = base_spec(outputs=[output])
            spec["source"]["fit_class"] = "both"
            self.assertInvalid(spec, "one fit per sample")

    def test_outputs_required(self):
        self.assertInvalid(base_spec(outputs=[]), "at least one output")

    def test_patient_ids_rejected(self):
        self.assertInvalid(base_spec(samples={"source": "list", "ids": ["-".join(["P", "1".zfill(7)])]}), "Patient IDs")

    def test_empty_sample_list_rejected(self):
        self.assertInvalid(base_spec(samples={"source": "list", "ids": []}), "at least one sample")

    def test_range_filters(self):
        self.assertInvalid(base_spec(filters={"purity": {}}), "minimum or a maximum")
        self.assertInvalid(base_spec(filters={"purity": {"min": "0.3"}}), "must be a number")
        self.assertInvalid(base_spec(filters={"purity": {"min": True}}), "must be a number")
        self.assertInvalid(base_spec(filters={"purity": {"min": float("nan")}}), "must be a number")
        self.assertInvalid(base_spec(filters={"purity": {"min": 30}}), "between")
        self.assertInvalid(base_spec(filters={"ploidy": {"min": 5, "max": 1}}), "above the maximum")

    def test_clinical_filter_needs_clinical_file(self):
        spec = base_spec(filters={"tmb": {"min": 1}})
        spec["source"]["clinical_file"] = ""
        self.assertInvalid(spec, "need a clinical file")

    def test_clinical_purity_fraction_warns(self):
        _, warnings = S.validate(base_spec(filters={"clinical_purity": {"max": 1}}))
        self.assertTrue(any("percentage" in w for w in warnings))

    def test_scope_warnings(self):
        _, warnings = S.validate(base_spec(fits={"single_run": False, "tiers": ["default"]}))
        self.assertTrue(any("every fit" in w for w in warnings))
        _, warnings = S.validate(base_spec(outputs=[{"type": "gene_level"}]))
        self.assertTrue(any("all genes" in w for w in warnings))
        _, warnings = S.validate(base_spec(outputs=[{"type": "copy_folders"}]))
        self.assertTrue(any("disk space" in w for w in warnings))

    def test_suggest_resources(self):
        self.assertEqual(S.suggest_resources(base_spec(samples={"source": "list", "ids": ["a"]})), "small")
        self.assertEqual(S.suggest_resources(base_spec(samples={"source": "list", "ids": ["a"] * 501})), "medium")
        self.assertEqual(S.suggest_resources(base_spec()), "medium")
        self.assertEqual(S.suggest_resources(base_spec(fits={"single_run": False, "tiers": ["default"]})), "large")
        self.assertEqual(S.suggest_resources(base_spec(outputs=[{"type": "gene_level"}])), "large")
        self.assertEqual(S.suggest_resources(base_spec(outputs=[{"type": "copy_folders"}])), "large")


class RenderTests(unittest.TestCase):
    def test_invalid_spec_refused(self):
        with self.assertRaises(ValueError):
            S.render(base_spec(outputs=[]), ctx())

    def test_every_option_compiles(self):
        compile(S.render(every_option_spec(), ctx()), "script.py", "exec")

    def test_each_output_alone_compiles(self):
        for output in every_option_spec()["outputs"]:
            compile(S.render(base_spec(outputs=[output]), ctx()), "script.py", "exec")

    def test_call_order(self):
        script = S.render(every_option_spec(), ctx())
        order = ["setFitPolicy", "selectSamplesFromFile", "meta.buildFacetsMeta", "setPurityFilter",
                 "dataset.buildFacetsDataset", "runAlterationAnalysis", "writeAlterationData"]
        positions = [script.index(call) for call in order]
        self.assertEqual(positions, sorted(positions))

    def test_human_fits_only(self):
        script = S.render(base_spec(fits={"single_run": True, "tiers": ["reviewed"]}), ctx())
        self.assertIn("meta.setSingleRunPerSample(True, False)", script)
        self.assertIn("meta.setFitPolicy(FitPolicy.reviewed_only())", script)

    def test_segment_scores_run_before_alteration_analysis(self):
        script = S.render(every_option_spec(), ctx())
        for call in ("cin_for_runs(", "chromothripsis_for_runs(", "compute_sample_percentages("):
            self.assertLess(script.index(call), script.index("runAlterationAnalysis"), call)

    def test_tie_break(self):
        script = S.render(base_spec(fits={"single_run": True, "tiers": ["reviewed"], "tie_break": "newest_reviewed"}), ctx())
        self.assertIn("FitPolicy.reviewed_only(tie_break='newest_reviewed')", script)
        script = S.render(base_spec(fits={"single_run": True, "tiers": ["reviewed", "autoqc"], "tie_break": "newest_reviewed"}), ctx())
        self.assertIn("FitPolicy(('reviewed', 'autoqc'), tie_break='newest_reviewed')", script)

    def test_ascets_uses_the_servers_r_prefix(self):
        script = S.render(base_spec(outputs=[{"type": "ascets"}]), ctx(r_prefix="module load r/3.6.3"))
        self.assertIn("FacetsConfig.from_env(r_prefix='module load r/3.6.3')", script)

    def test_numbers_stay_numbers(self):
        script = S.render(base_spec(filters={"purity": {"min": 0.3}, "tmb": {"min": 10}}), ctx())
        self.assertIn("dataset.setPurityFilter(0.3, 1.0)", script)
        self.assertIn('dataset.setTMBFilter(10, float("inf"))', script)

    def test_injection_stays_inert(self):
        marker = os.path.join(tempfile.gettempdir(), "facetsapi_spec_injected")
        payload = '"); open(%r, "w"); ("' % marker
        spec = base_spec(filters={"cancer_type": [payload]})
        spec["source"]["facets_dir"] = "/repo/" + payload
        script = S.render(spec, ctx(output_dir="/out/" + payload))
        namespace_calls = [line for line in script.splitlines() if line.startswith("open(")]
        self.assertEqual(namespace_calls, [])
        compile(script, "script.py", "exec")

    def test_newlines_rejected_or_neutralised(self):
        code = "x\nimport shutil"
        self.assertTrue(S.validate(base_spec(name=code))[0])
        self.assertTrue(S.validate(base_spec(filters={"cancer_type": [code]}))[0])
        spec = base_spec()
        spec["source"]["facets_dir"] = "/repo/" + code
        self.assertTrue(S.validate(spec)[0])
        self.assertTrue(S.validate(base_spec(samples={"source": "list", "ids": [code]}))[0])
        self.assertTrue(S.validate(base_spec(outputs=[{"type": "gene_level", "genes": [code]}]))[0])
        script = S.render(base_spec(), ctx(user="me\nimport shutil", job_id="1\nimport shutil", facetsapi_sha="a\nimport shutil"))
        self.assertNotIn("\nimport shutil", script)

    def test_samples_text(self):
        self.assertIsNone(S.samples_text(base_spec()))
        self.assertEqual(S.samples_text(base_spec(samples={"source": "list", "ids": ["a", "b"]})), "a\nb\n")


class RepositoryCase(unittest.TestCase):
    """Render scripts and run them in a fresh interpreter against the synthetic repository."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.info = build(cls.tmp.name)
        cls.run_count = 0

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def fake_rscript(self):
        """Stands in for R: run_ascets.R needs tidyverse, which is not a test dependency."""
        path = os.path.join(self.tmp.name, "fake_rscript")
        with open(path, "w") as fh:
            fh.write('#!/bin/bash\necho "$@" > "$(dirname "$0")/rscript_args"\n')
        os.chmod(path, 0o755)
        return path

    def repo_snapshot(self):
        """Every path in the repository with its mode, size, mtime and content."""
        snapshot = {}
        for root, dirs, files in os.walk(self.info["repo"]):
            for name in dirs + files:
                path = os.path.join(root, name)
                info = os.lstat(path)
                digest = None
                if os.path.isfile(path):
                    with open(path, "rb") as fh:
                        digest = hashlib.sha256(fh.read()).hexdigest()
                snapshot[path] = (info.st_mode, info.st_size, info.st_mtime_ns, digest)
        return snapshot

    def run_spec(self, spec, samples=None, edit_script=None, expect_failure=False, real_r=False):
        type(self).run_count += 1
        job_dir = os.path.join(self.tmp.name, "job%d" % self.run_count)
        out_dir = os.path.join(job_dir, "out")
        self.last_job_dir = job_dir
        os.makedirs(job_dir)
        spec = copy.deepcopy(spec)
        spec["source"]["facets_dir"] = self.info["repo"]
        spec["source"]["clinical_file"] = self.info["clinical_file"]
        if samples:
            spec["samples"] = {"source": "list", "ids": samples}
        text = S.samples_text(spec)
        if text:
            with open(os.path.join(job_dir, "samples.txt"), "w") as fh:
                fh.write(text)
        script = os.path.join(job_dir, "script.py")
        with open(script, "w") as fh:
            rendered = S.render(spec, ctx(output_dir=out_dir, job_dir=job_dir, protected_roots=[self.info["repo"]]))
            fh.write(edit_script(rendered) if edit_script else rendered)
        env = dict(os.environ) if real_r else dict(os.environ, FACETSAPI_RSCRIPT=self.fake_rscript())
        result = subprocess.run([sys.executable, script], capture_output=True, text=True, cwd=job_dir, env=env)
        if expect_failure:
            self.assertNotEqual(result.returncode, 0, result.stdout[-2000:])
            return result
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr[-2000:])
        with open(os.path.join(job_dir, "run_summary.json")) as fh:
            return out_dir, json.load(fh)

    def rows(self, path):
        with open(path) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


class EndToEndTests(RepositoryCase):
    def test_report_and_gene_level(self):
        tumor = self.info["standard"]["reviewed_best"][0]
        out, summary = self.run_spec(base_spec(outputs=[{"type": "report"}, {"type": "gene_level", "genes": ["TP53"]}]),
                                     samples=[tumor])
        self.assertEqual(summary["num_runs"], 1)
        self.assertEqual(len(self.rows(os.path.join(out, "report.txt"))), 1)
        self.assertTrue(os.path.exists(os.path.join(out, "gene_level_report.txt")))
        self.assertEqual(sorted(summary["outputs"]), ["gene_level_report.txt", "report.txt"])

    def test_reviewed_only_drops_samples_without_review(self):
        samples = [self.info["standard"]["reviewed_best"][0], self.info["standard"]["autoqc_only"][0]]
        _, everything = self.run_spec(base_spec(), samples=samples)
        _, reviewed = self.run_spec(base_spec(fits={"single_run": True, "tiers": ["reviewed"]}), samples=samples)
        self.assertEqual(everything["num_runs"], 2)
        self.assertEqual(reviewed["num_runs"], 1)

    def test_unmatched_ids_reported(self):
        tumor = self.info["standard"]["reviewed_best"][0]
        _, summary = self.run_spec(base_spec(), samples=[tumor, sid(9999999)])
        self.assertEqual(summary["unmatched_sample_ids"], [sid(9999999)])

    def test_filter_with_one_sided_bound(self):
        _, kept = self.run_spec(base_spec(filters={"tmb": {"min": 0}}))
        _, dropped = self.run_spec(base_spec(filters={"tmb": {"min": 100000}}))
        self.assertGreater(kept["num_runs"], 0)
        self.assertEqual(dropped["num_runs"], 0)

    def test_clinical_purity_uses_percent(self):
        _, summary = self.run_spec(base_spec(filters={"clinical_purity": {"min": 0}}))
        self.assertGreater(summary["num_runs"], 0)

    def test_every_output(self):
        out, summary = self.run_spec(base_spec(outputs=every_option_spec()["outputs"]),
                                     samples=[self.info["standard"]["reviewed_best"][0], self.info["standard"]["autoqc_only"][0]])
        self.assertGreater(summary["num_runs"], 0)
        for name in ("report.txt", "gene_level_report.txt", "alteration_data.txt", "purity_cf.txt",
                     "merged_adjusted.seg", "merged_cncf.txt", "cin_full.tsv", "cin_arm.tsv", "cin_chrom.tsv",
                     "chromothripsis_cin_full.tsv", "chromothripsis_cin_arm.tsv", "chromothripsis_candidate_arms.tsv",
                     "hetloss_homdel_percent.tsv", "ascets/merged_seg_data.seg", "facets_folders/all"):
            self.assertTrue(os.path.exists(os.path.join(out, name)), name)
        self.assertEqual(len(glob.glob(os.path.join(out, "FacetsDataset_summary_*.tsv"))), 1)
        self.assertEqual(len(glob.glob(os.path.join(out, "FacetsDataset_histogram_*.png"))), 2)

    def test_scores_match_the_standalone_tools(self):
        samples = [self.info["standard"]["reviewed_best"][0], self.info["standard"]["autoqc_only"][0]]
        out, _ = self.run_spec(base_spec(outputs=[{"type": "cin", "levels": ["full", "arm"]}, {"type": "alteration"}]), samples=samples)
        code = ("import sys\nsys.path[:0] = [%r, %r]\nfrom cin_scores import calculate_cin, write_cin_results\n"
                "for level, split in (('full', True), ('arm', True)):\n"
                "    write_cin_results(calculate_cin(%r, method='Eckert', splitArms=split, cinLevel=level, "
                "clinical_sample_file=%r, facets_dir=%r), %r + '/' + level + '.tsv')\n"
                % (CHECKOUT, os.path.join(CHECKOUT, "scripts"), samples, self.info["clinical_file"], self.info["repo"], out))
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr[-1500:])
        for level in ("full", "arm"):
            self.assertEqual(open(os.path.join(out, "cin_%s.tsv" % level)).read(), open(os.path.join(out, level + ".tsv")).read())

    def test_cin_leaves_runs_unchanged_and_skips_missing_ploidy(self):
        code = ("import sys\nsys.path[:0] = [%r, %r]\nfrom facetsAPI import FacetsMeta, FacetsDataset\n"
                "from cin_scores import cin_for_runs\n"
                "meta = FacetsMeta(%r, %r, 'purity'); meta.setSingleRunPerSample(True, True)\n"
                "meta.build_from_file_listing = True; meta.samples_from_file.append(%r); meta.buildFacetsMeta()\n"
                "ds = FacetsDataset(meta); ds.buildFacetsDataset(); run = ds.runList[0]\n"
                "before = len(run.segments); arm = cin_for_runs([run], cinLevel='arm')\n"
                "assert len(run.segments) == before and len(arm) == 1\n"
                "run.ploidy = 'NA'; assert cin_for_runs([run]) == []\n"
                % (CHECKOUT, os.path.join(CHECKOUT, "scripts"), self.info["clinical_file"], self.info["repo"],
                   self.info["standard"]["reviewed_best"][0]))
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr[-1500:])
        self.assertIn("skipped 1 runs", result.stdout)

    def test_ascets_reaches_r_with_chosen_options(self):
        tumor = self.info["standard"]["reviewed_best"][0]
        out, _ = self.run_spec(base_spec(outputs=[{"type": "ascets", "genome": "hg38", "min_arm_breadth": 0.4,
                                                   "arm_alt_frac_thresh": 0.6, "keep_noise": True, "use_adjusted": False}]),
                               samples=[tumor])
        args = open(os.path.join(self.tmp.name, "rscript_args")).read()
        self.assertIn("run_ascets.R", args)
        self.assertIn(os.path.join(out, "ascets", ""), args)  # the R code appends file names with no separator
        self.assertIn("genomic_arm_coordinates_hg38.txt", args)
        for piece in ("-m 0.4", "-k T", "-a 0.6", "facetsUnadjustedSeg"):
            self.assertIn(piece, args)
        self.assertTrue(os.path.exists(os.path.join(out, "ascets", "merged_seg_data.seg")))

    def test_merged_files_and_ascets_respect_filters(self):
        samples = [self.info["standard"]["reviewed_best"][0], self.info["standard"]["autoqc_only"][0]]
        outputs = [{"type": "report"}, {"type": "merged", "files": ["adjusted_seg", "cncf"]}, {"type": "ascets"}]
        out, _ = self.run_spec(base_spec(outputs=outputs), samples=samples)
        unfiltered = open(os.path.join(out, "merged_adjusted.seg")).read()
        self.assertIn(samples[0], unfiltered)
        self.assertIn(samples[1], unfiltered)

        tmb = self.rows(os.path.join(out, "report.txt"))
        keep = max(tmb, key=lambda row: float(row["tmb"]))
        drop = min(tmb, key=lambda row: float(row["tmb"]))
        self.assertNotEqual(float(keep["tmb"]), float(drop["tmb"]))
        spec = base_spec(outputs=outputs, filters={"tmb": {"min": float(keep["tmb"])}})
        out, summary = self.run_spec(spec, samples=samples)
        self.assertEqual(summary["num_runs"], 1)
        for path in ("merged_adjusted.seg", "merged_cncf.txt", os.path.join("ascets", "merged_seg_data.seg")):
            text = open(os.path.join(out, path)).read()
            self.assertIn(keep["ID"], text, path)
            self.assertNotIn(drop["ID"], text, path)

    def test_ascets_with_no_samples_says_so(self):
        spec = base_spec(outputs=[{"type": "ascets"}], filters={"tmb": {"min": 100000}})
        result = self.run_spec(spec, expect_failure=True)
        self.assertIn("no samples passed", result.stderr)

    def test_all_fits_mode(self):
        _, single = self.run_spec(base_spec())
        _, every = self.run_spec(base_spec(fits={"single_run": False, "tiers": ["default"]}))
        self.assertGreater(every["num_runs"], single["num_runs"])

    def test_2n_fit_classes(self):
        tumor = self.info["twon"]["full"][0]
        for fit_class, expected in (("clinical", 1), ("research", 1), ("both", 2)):
            spec = base_spec()
            spec["source"]["fit_class"] = fit_class
            _, summary = self.run_spec(spec, samples=[tumor])
            self.assertEqual(summary["num_runs"], expected, fit_class)

    def test_injection_payload_runs_inert(self):
        marker = os.path.join(self.tmp.name, "injected")
        payload = "'); open(%r, 'w'); ('" % marker
        out, summary = self.run_spec(base_spec(filters={"cancer_type": [payload]}))
        self.assertEqual(summary["num_runs"], 0)
        self.assertFalse(os.path.exists(marker))


class RepositoryIsNeverChangedTests(RepositoryCase):
    """The hard rule: nothing a job does may write, modify or delete anything in a repository."""

    def assert_repo_untouched_by(self, before):
        self.assertEqual(self.repo_snapshot(), before)

    def test_every_output_leaves_the_repository_identical(self):
        before = self.repo_snapshot()
        self.assertGreater(len(before), 50)
        out, summary = self.run_spec(base_spec(outputs=every_option_spec()["outputs"]))
        self.assertGreater(summary["num_runs"], 0)
        self.assertTrue(os.path.exists(os.path.join(out, "facets_folders", "all")))
        self.assert_repo_untouched_by(before)

    def blocked_when_added_to_script(self, code):
        before = self.repo_snapshot()
        result = self.run_spec(base_spec(), edit_script=lambda s: s.replace('dataset.writeReport', code + "\ndataset.writeReport", 1),
                               expect_failure=True)
        self.assertIn("Blocked", result.stderr, result.stderr[-1500:])
        self.assert_repo_untouched_by(before)

    def test_backup_methods_are_blocked(self):
        self.blocked_when_added_to_script("dataset.backupMafs(meta)")

    def test_writing_a_repository_file_is_blocked(self):
        self.blocked_when_added_to_script("open(os.path.join(meta.facets_repo_path, 'x.txt'), 'w')")

    def test_appending_to_a_repository_file_is_blocked(self):
        self.blocked_when_added_to_script("open(dataset.runList[0].fitDir + '/../facets_qc.txt', 'a')")

    def test_deleting_is_blocked(self):
        self.blocked_when_added_to_script("os.remove(os.path.join(dataset.runList[0].fitDir, os.listdir(dataset.runList[0].fitDir)[0]))")

    def test_deleting_a_tree_is_blocked(self):
        self.blocked_when_added_to_script("import shutil; shutil.rmtree(dataset.runList[0].fitDir)")

    def test_renaming_out_of_a_repository_is_blocked(self):
        self.blocked_when_added_to_script("os.rename(dataset.runList[0].fitDir, os.path.join(OUT, 'stolen'))")

    def test_chmod_is_blocked(self):
        self.blocked_when_added_to_script("os.chmod(dataset.runList[0].fitDir, 0o777)")

    def test_symlink_from_output_into_repository_is_blocked(self):
        self.blocked_when_added_to_script("os.symlink(meta.facets_repo_path, os.path.join(OUT, 'link')); "
                                          "open(os.path.join(OUT, 'link', 'x.txt'), 'w')")

    def test_starting_a_program_is_blocked(self):
        self.blocked_when_added_to_script("import subprocess; subprocess.run(['touch', os.path.join(meta.facets_repo_path, 'x')])")

    def test_shell_is_blocked(self):
        self.blocked_when_added_to_script("os.system('touch ' + os.path.join(meta.facets_repo_path, 'x'))")

    def test_a_faked_ascets_command_needs_ascets_selected(self):
        self.blocked_when_added_to_script("import subprocess; subprocess.run(['bash', '-lc', 'rm -rf x run_ascets.R'])")

    @unittest.skipUnless(r_has_ascets_packages(), "needs Rscript with dplyr, tidyr, tibble and data.table")
    def test_real_ascets_writes_only_to_its_output_folder(self):
        before = self.repo_snapshot()
        samples = [self.info["standard"]["reviewed_best"][0], self.info["standard"]["autoqc_only"][0]]
        out, _ = self.run_spec(base_spec(outputs=[{"type": "ascets", "min_arm_breadth": 0.1}]), samples=samples, real_r=True)
        names = os.listdir(os.path.join(out, "ascets"))
        for suffix in ("_arm_level_calls.txt", "_aneuploidy_scores.txt", "_arm_weighted_average_segmeans.txt"):
            self.assertTrue(any(n.endswith(suffix) for n in names), (suffix, names))
        calls = next(n for n in names if n.endswith("_arm_level_calls.txt"))
        text = open(os.path.join(out, "ascets", calls)).read()
        self.assertIn(samples[0], text)
        self.assertEqual(sorted(os.listdir(out)), ["ascets"])
        self.assert_repo_untouched_by(before)

    def test_output_folder_symlinked_into_repository_refuses_to_start(self):
        before = self.repo_snapshot()
        link = os.path.join(self.tmp.name, "innocent_looking")
        os.symlink(self.info["repo"], link)
        job_dir = os.path.join(self.tmp.name, "symlink_job")
        os.makedirs(job_dir)
        script = os.path.join(job_dir, "script.py")
        spec = base_spec()
        spec["source"]["facets_dir"] = self.info["repo"]
        spec["source"]["clinical_file"] = self.info["clinical_file"]
        with open(script, "w") as fh:
            fh.write(S.render(spec, ctx(output_dir=os.path.join(link, "out"), job_dir=job_dir, protected_roots=[self.info["repo"]])))
        result = subprocess.run([sys.executable, script], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("overlaps a protected repository", result.stderr)
        self.assert_repo_untouched_by(before)


class CliTests(unittest.TestCase):
    def run_cli(self, spec, command):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "spec.json")
            with open(path, "w") as fh:
                json.dump(spec, fh)
            return subprocess.run([sys.executable, "-m", "facetsAPI.spec", command, path], capture_output=True,
                                  text=True, cwd=CHECKOUT)

    def test_render(self):
        result = self.run_cli(base_spec(), "render")
        self.assertEqual(result.returncode, 0, result.stderr)
        compile(result.stdout, "script.py", "exec")

    def test_validate_reports_errors(self):
        result = self.run_cli(base_spec(outputs=[]), "validate")
        self.assertEqual(result.returncode, 1)
        self.assertIn("at least one output", result.stderr)


class StandaloneLoadTests(unittest.TestCase):
    def test_loads_by_path_without_pandas(self):
        code = ("import sys, importlib.util as u\n"
                "sys.modules['pandas'] = None\n"
                "s = u.spec_from_file_location('fspec', %r)\n"
                "m = u.module_from_spec(s); s.loader.exec_module(m)\n"
                "print(m.SPEC_VERSION)\n" % os.path.join(CHECKOUT, "facetsAPI", "spec.py"))
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()

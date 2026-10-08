from heapq import merge
from random import sample
import sys
import os.path
import os
import glob
import shutil
import subprocess
import shlex
import pandas as pd
import numpy as np
import pickle
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
import statistics
from datetime import datetime
from datetime import date
import traceback


def _pyplot():
    """matplotlib is only needed by FacetsDataset.createHistogram(); import it lazily so the package
    imports on systems without it (e.g. a bare cluster python)."""
    try:
        import matplotlib
    except ImportError as e:
        raise FacetsConfigError("createHistogram() needs matplotlib: pip install matplotlib") from e
    if not os.environ.get("MPLBACKEND") and not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
        matplotlib.use("Agg")  # headless (cluster / container) default; MPLBACKEND still wins
    import matplotlib.pyplot as plt
    return plt

from .exceptions import (FacetsAPIError, FacetsConfigError, FacetsRepositoryError,
                         ManifestError, ParseError, ExternalToolError)
from .cn_calls import CN_CALL_TABLE as _CN_CALL_TABLE, get_cn_call as _get_cn_call
from .layout import (discover_sample, select_classes, read_manifest, short_sample_id, pair_id, is_ultra,
                     fit_dir as _fit_dir, locate_fit_files, validate_fit_class, validate_mode, NOT_SELECTED)
from .fitpolicy import (FitPolicy, FitCandidate, tier_of, label_of,
                        TIER_REVIEWED, TIER_AUTOQC, TIER_DEFAULT)
from .config import FacetsConfig


def _to_float(value):
    """None / NaN / 'NA' / '' -> None, otherwise float(value) or None when unparseable."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, str) and value.strip().upper() in ("", "NA", "NAN", "NONE"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_bool(value):
    """Interpret facets TRUE/FALSE strings, numpy bools and plain bools; None/NaN -> False."""
    if isinstance(value, str):
        return value.strip().upper() in ("TRUE", "T", "1", "YES")
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    return bool(value)



######################
# bcolors:    This class is a set of color codes that can be concatenated in strings to change output in the terminal.
#             Can start with the color and must end with the ENDC color.
######################
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


######################
# ExtTools:    This class functions as a wrapper to execute certain functions from external sources.
#              R is invoked through FacetsConfig (see config.py): set FACETSAPI_R_PREFIX / FACETSAPI_RSCRIPT
#              or pass config=FacetsConfig(...) to control how Rscript is found. Nothing is submitted to a
#              scheduler; wrap the call yourself (sbatch, a Nextflow process, a container) if you need one.
######################
class ExtTools:
    def __init__(self, ref_meta, config=None):
        if not isinstance(ref_meta, FacetsMeta):
            raise FacetsConfigError("ExtTools(): ref_meta must be of type FacetsMeta.")
        self.ref_meta = ref_meta
        self.config = config if config is not None else FacetsConfig.from_env()
        if not isinstance(self.config, FacetsConfig):
            raise FacetsConfigError("ExtTools(): config must be a FacetsConfig.")

    #Running FACETS itself is the job of the processing pipeline (facets-suite run-facets-wrapper.R), not this API.
    def runFacets(self):
        raise NotImplementedError("ExtTools.runFacets() is not implemented; FACETS runs are produced by the IMPACT "
                                  "FACETS pipeline (facets-suite run-facets-wrapper.R), which this API only reads.")

    #This function will make a merged file for any file type in the MetaDictMap with a shared header.
    #fileType should be a MetaDictMap.FILETYPE.
    #This will return a path to the written file.
    def makeMergedFile(self, fileType, outfile_path):
        #Keep track and don't duplicate any files.
        #This can happen when there are alt fits included for the same sample.
        #For instance a seg file merge could merge a default and alt fit, but would duplicate
        #the exact values.
        completed_files = set()
        wrote_header = False
        with open(outfile_path, 'w') as outfile:
            for item in self.ref_meta.master_file_dict:
                cur_file = self.ref_meta.master_file_dict.get(item)[fileType]
                if not cur_file:
                    continue
                base_file = os.path.basename(cur_file)
                if base_file in completed_files:
                    continue
                completed_files.add(base_file)
                with open(cur_file) as infile:
                    header = infile.readline()
                    if not wrote_header:
                        outfile.write(header if header.endswith("\n") else header + "\n")
                        wrote_header = True
                    for line in infile:
                        outfile.write(line if line.endswith("\n") else line + "\n")
        return outfile_path

    #This function will run ascets (arm-level aneuploidy calls) over the merged seg file of this object's FacetsMeta.
    #It runs inline (no scheduler) and writes runLog.out / runLog.err into output_dir. Returns the output prefix.
    def runAscets(self, output_dir="", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True):
        if output_dir == "":
            raise FacetsConfigError("ExtTools.runAscets(): No output path provided.")
        if ref_genome_coords not in ("hg19", "hg38"):
            raise FacetsConfigError("ExtTools.runAscets(): ref_genome_coords must be hg19 or hg38.")
        arm_coords_file   = os.path.join(self.config.ascets_dir, "genomic_arm_coordinates_%s.txt" % ref_genome_coords)
        run_ascets_script = os.path.join(self.config.ascets_dir, "run_ascets.R")
        for needed in (arm_coords_file, run_ascets_script):
            if not os.path.isfile(needed):
                raise FacetsConfigError("ExtTools.runAscets(): missing %s (set FACETSAPI_ASCETS_DIR or FacetsConfig.ascets_dir)." % needed)

        #Print some notes on dependencies.
        print("\t\tNote: Ascets requires R 3.6.* and requires tidyverse and data.table libraries.")
        print("\t\tNote: Failure to meet these dependencies will result in an empty ascets run.")

        #Make an output folder if necessary.
        if os.path.isdir(output_dir):
            print("\t\tAscets folder already exists: " + output_dir)
        else:
            print("\t\tCreating ascets output folder at: " + output_dir)
            os.makedirs(output_dir, exist_ok=True)

        #Make an adjusted or unadjusted merged seg file to use as ascets input.
        #This will use everything in the FacetsMeta reference object for ExtTools.
        merged_seg_tmp_file = os.path.join(output_dir, "merged_seg_data.seg")
        if use_adjusted:
            adj_str = "facetsAdjustedSeg"
            self.makeMergedFile(MetaDictMap.ADJUSTED_SEG_FILE, merged_seg_tmp_file)
        else:
            adj_str = "facetsUnadjustedSeg"
            self.makeMergedFile(MetaDictMap.UNADJUSTED_SEG_FILE, merged_seg_tmp_file)

        keep_noise_short = "T" if keep_noise else "F"
        noise_str        = "keepNoise" if keep_noise else "filterNoise"

        #Make an output prefix:
        out_prefix = "FacetsAscetsRun_" + ref_genome_coords + "_breadth" + str(min_arm_breadth) + "_" + noise_str + "_armAltFracThresh" + str(arm_alt_frac_thresh) + "_" + adj_str + "_" + str(date.today())

        args = [run_ascets_script, "-i", merged_seg_tmp_file, "-c", arm_coords_file, "-m", str(min_arm_breadth),
                "-k", keep_noise_short, "-a", str(arm_alt_frac_thresh), "-o", out_prefix, "-p", output_dir]
        print("\t\t\tRunning ascets...")
        print("\t\t\t" + self.config.compose_r_command(args))
        result = self.config.run_r(args, check=False)
        with open(os.path.join(output_dir, "runLog.out"), "w") as fh:
            fh.write(result.stdout or "")
        with open(os.path.join(output_dir, "runLog.err"), "w") as fh:
            fh.write(result.stderr or "")
        if result.returncode != 0:
            raise ExternalToolError(self.config.compose_r_command(args), result.returncode, result.stderr)
        return out_prefix


######################
# FPTools:    This class functions as a wrapper to execute certain functions from Facets Preview.
#             The facets-preview config json is required (config_file=..., FACETSAPI_FP_CONFIG, or FacetsConfig.fp_config_file).
######################
class FPTools:
    def __init__(self, config_file=None, cbio_maf_file="", cbio_nonsigned_maf_file="", config=None):
        self.config = config if config is not None else FacetsConfig.from_env()
        if not isinstance(self.config, FacetsConfig):
            raise FacetsConfigError("FPTools(): config must be a FacetsConfig.")
        if config_file:
            self.config = self.config.with_(fp_config_file=config_file)
        self.cbio_maf_file           = cbio_maf_file
        self.cbio_nonsigned_maf_file = cbio_nonsigned_maf_file

    @property
    def config_file(self):
        return self.config.fp_config_file

    @property
    def genomic_anno_wrapper_script(self):
        return self.config.fp_wrapper_script

    #This function runs the wrapper script, R/fp_generate_genomic_anno_wrapper.R
    #That script takes these arguments in this order:
    #facets_preview_config_file = args[1]
    #tumor_sample = args[2] -- i.e. P-0000280-T05-IM7
    #sample_id = args[3] -- i.e. P-0071250-T01-IM7_P-0071250-N01-IM7
    #sample_path = args[4]
    #cbio_maf_path = args[5]
    #cbio_maf_nonsignedout_path = args[6]
    #r_lib_loc = args[7] (optional; the R library holding facetsPreview)
    def runGenerateGenomicAnnotations(self, short_id, long_id, path):
        if not self.config.fp_config_file:
            raise FacetsConfigError("FPTools.runGenerateGenomicAnnotations(): facets-preview config file not set. "
                                    "Pass FPTools(config_file=...) or set FACETSAPI_FP_CONFIG.")
        if not os.path.isfile(self.config.fp_wrapper_script):
            raise FacetsConfigError("FPTools.runGenerateGenomicAnnotations(): wrapper script not found: %s "
                                    "(set FACETSAPI_FP_WRAPPER)." % self.config.fp_wrapper_script)
        args = [self.config.fp_wrapper_script, self.config.fp_config_file, short_id, long_id, path,
                self.cbio_maf_file, self.cbio_nonsigned_maf_file]
        if self.config.r_lib_loc:
            args.append(self.config.r_lib_loc)
        print(self.config.compose_r_command(args))
        return self.config.run_r(args)

    @staticmethod
    def _backupFiles(ref_meta, slot, what):
        if not isinstance(ref_meta, FacetsMeta):
            raise FacetsConfigError("FPTools backup: ref_meta must be of type FacetsMeta.")
        print("Backing up " + what + " files.")
        done = set()
        for item in ref_meta.master_file_dict:
            cur_file = ref_meta.master_file_dict.get(item)[slot]
            if not cur_file or cur_file in done or not os.path.isfile(cur_file):
                continue
            done.add(cur_file)
            shutil.copy2(cur_file, cur_file + "." + str(date.today()) + ".bak")
        return len(done)

    #This function will make backup files of existing manifests for every sample in the provided FacetsMeta object.
    def backupManifests(self, ref_meta):
        return FPTools._backupFiles(ref_meta, MetaDictMap.MANIFEST_FILE, "manifest")

    #This function will make backup files of existing facets qc files for every sample in the provided FacetsMeta object.
    def backupFacetsQC(self, ref_meta):
        return FPTools._backupFiles(ref_meta, MetaDictMap.FACETS_QC_FILE, "facets_qc")


class MetaDictMap:
    OUT_FILE = 0
    CNCF_FILE = 1
    QC_FILE = 2
    FACETS_QC_FILE = 3
    SELECTED_FIT_DIR = 4
    GENE_FILE = 5
    ADJUSTED_SEG_FILE = 6
    CCF_MAF_FILE = 7
    NONSIGNED_MAF_FILE = 8
    UNADJUSTED_SEG_FILE = 9
    MANIFEST_FILE = 10
    FIT_STATUS = 11
    SAMPLE_BASE_DIR = 12   # the discovered sample dir (standard) or class dir (2n) that holds the manifest
    RDATA_FILE = 13        # "" when the fit has no .Rdata
    FIT_CLASS = 14         # "standard" | "clinical" | "research"
    FIT_TIER = 15          # "reviewed" | "autoqc" | "default" | "unreviewed" (tier the fit was selected/loaded under)

######################
# FacetsMeta:    This class represents necessary metadata structures that hold information such as directory paths,
#                cancer types, and other clinical data that need to be referenced in FacetsRuns.
######################
class FacetsMeta:
    use_unreviewed_defaults = False
    hisens_vs_purity        = "purity"
    selectSingleRun         = False
    failed_samples          = []
    run_verbose             = False
    build_from_file_listing = False

    #OncoTree Codes that correspond to a specific cancer type.  
    breast_carcinoma    = ["ILC","IDC","BRCA","BRCNOS","BRCANOS","MDLC","MBC","CSNOS"]
    nscLung_cancer      = ["LUAD","LUSC","LUNE","NSCLC","LUAS","NSCLCPD","ALUCA","SARCL"]
    endometrial_cancer  = ["UMEC","UCCC","UEC","UCEC","UCS","UDDC","UUC","USC","OUTT"]
    colorectal_cancer   = ["READ","COAD","MACR","COADREAD","CMC","CM"]
    melanoma            = ["SKCM","ACRM","MUP","ARMM","HNMUCM","UM","VMM","SKCN"]
    scLung_cancer       = ["SCLC","CSCLC"]
    ovarian_cancer      = ["HGSOC","CCOV","LGSOC","OCS","EOV"]
    prostate_cancer     = ["PRAD","PRSCC","PRSC"]
    unknown_primary     = ["NETNOS","ADNOS","CUP","PDC","NECNOS","CUPNOS","SCUP"]
    bladder_cancer      = ["BLCA", "UTUC", "UCU", "BLSC", "USCC"]
    pancreatic_cancer   = ["PAAD","PANET","PAASC","UCP","SPN","PB"]
    renalCell_carcinoma = ["PRCC","TRCC","CCRCC","URCC","RCC","CHRCC","MT","SRCC","ROCY","MTSCC"]
    glioma              = ["ASTR","ODG","AODG","GBM","HGGNOS","AASTR","GB","DIFG"]
    headNeck_carcinoma  = ["HNSC","OPHSC","HPHSC","OCSC","HNSCUP","LXSC","HNNE","SNSC","ODGC"]
    germCell_tumor      = ["MGCT","SEM","VMT","BMT","OYST","GCTSTM","NSGCT","EMBCA","OGCT","OMGCT","TT","TYST","VDYS","ODYS","OIMT","VYST","VMGCT","BMGCT","BIMT","VIMT","BYST","OMT","GCT"]

    # Copy-number call table keyed on (wgd, tcn, lcn). Lives in cn_calls.py; aliased here
    # so existing code referring to FacetsMeta.CN_CALL_TABLE keeps working.
    CN_CALL_TABLE = _CN_CALL_TABLE

    # Mapping of chromosome arms based on position in hg19. 
    # Format is chromosome: [p_start, p_end, q_start, q_end]
    # Information extracted in R using "data(hg19.chromArm, package="biscuiteer")".
    chr_arms = {
        1: [1,125000000,125000001,249250621],
        2: [1,93300000,93300001,243199373],
        3: [1,91000000,91000001,198022430],
        4: [1,50400000,50400001,191154276],
        5: [1,48400000,48400001,180915260],
        6: [1,61000000,61000001,171115067],
        7: [1,59900000,59900001,159138663],
        8: [1,45600000,45600001,146364022],
        9: [1,49000000,49000001,141213431],
        10: [1,40200000,40200001,135534747],
        11: [1,53700000,53700001,135006516],
        12: [1,35800000,35800001,133851895],
        13: [1,17900000,17900001,115169878],
        14: [1,17600000,17600001,107349540],
        15: [1,19000000,19000001,102531392],
        16: [1,36600000,36600001,90354753],
        17: [1,24000000,24000001,81195210],
        18: [1,17200000,17200001,78077248],
        19: [1,26500000,26500001,59128983],   
        20: [1,27500000,27500001,63025520],
        21: [1,13200000,13200001,48129895],
        22: [1,14700000,14700001,51304566],
        "X": [1,60600000,60600001,155270560],
        "Y": [1,12500000,12500001,59373566]
    }

    PERSIST_VERSION = 2   # bump when the pickled FacetsMeta state changes shape

    def __init__(self, clinical_sample_file = "", facets_repo_path = "", hisens_vs_purity="purity", persist_data="no",
                 fit_class="research", fit_policy=None):
        """
        clinical_sample_file : cbio data_clinical_sample(.oncokb).txt, or "" to walk the repository instead.
        facets_repo_path     : the FACETS repository root (the directory holding the P-00xxx folders).
        hisens_vs_purity     : which facets run to read per fit: "purity" (default) or "hisens".
        persist_data         : "no", or a path to pickle the built object to / load it from.
        fit_class            : for 2n repositories (pairs with clinical/ and research/ subtrees):
                               "clinical", "research" (default) or "both". Standard samples ignore it.
        fit_policy           : a FitPolicy controlling single-run fit selection. Default tiers:
                               reviewed > autoqc > default. See setSingleRunPerSample() / setFitPolicy().
        """
        #Relevant path and file data.
        self.clinical_sample_file    = clinical_sample_file if clinical_sample_file else ""
        self.facets_repo_path        = facets_repo_path if facets_repo_path else ""
        self.hisens_vs_purity        = validate_mode(hisens_vs_purity)
        self.persist_data            = persist_data if persist_data else "no"
        self.fit_class               = validate_fit_class(fit_class)
        if fit_policy is not None and not isinstance(fit_policy, FitPolicy):
            raise FacetsConfigError("fit_policy must be a FitPolicy instance.")
        self.fit_policy              = fit_policy if fit_policy is not None else FitPolicy()
        self._fit_policy_explicit    = fit_policy is not None

        #Per-instance run state. (Before 0.6.0 these were class attributes shared by every FacetsMeta.)
        self.selectSingleRun         = False
        self.use_unreviewed_defaults = False
        self.run_verbose             = False
        self.failed_samples          = []
        self.build_from_file_listing = False
        self.build_counts            = {}   # counters from the last parseClinicalSample()
        self._persist_version        = FacetsMeta.PERSIST_VERSION

        self.hasClinicalMeta = bool(self.clinical_sample_file)

        #Data structures and storage.

        self.master_file_dict       = {} # A map of relevant files for each sample. See MetaDictMap for ordering. {id -> [file list]}

        #These come from data_clinical_sample.
        self.cancer_type_map        = {} # A map of sample ids to cancer types.
        self.cancer_type_detail_map = {} # A map of sample ids to cancer type detailed.
        self.patient_id_map         = {} # A map of sample ids to patient ids.
        self.clinical_purity_map    = {} # A map of purity values from the clinical sample file.
        self.onkotree_code_map      = {} # A map of onkotree codes.
        self.cvr_tmb_score_map      = {} # A map of cvr tmb scores from the clinical sample file.
        self.msi_score_map          = {} # A map of msi scores.
        self.long_id_map            = {} # A map of sample ids to their corresponding long_ids.  id -> [long_id1, long_id2...]
        self.samples_from_file      = [] # A list of samples from a file that should be selected for this object.
        self.fit_map                = {} # A map of run key -> [fit label, fit dir]. Labels: Best, Acceptable, AutoQC_Best, AutoQC_Pass, Default, Unreviewed.

    @staticmethod
    def getArm(chrom, start, end):
            """
            Given a chromosome and a genomic interval [start, end],
            return the string '<chrom><arm>' where arm is 'p', 'q', or 'pq'
            (if the interval spans the centromere).
            """
            # normalize chromosome key for lookup in chr_arms
            try:
                key = int(chrom)
            except (ValueError, TypeError):
                key = chrom  # e.g. 'X' or 'Y'
            if key not in FacetsMeta.chr_arms:
                raise ValueError(f"getArm: unknown chromosome '{chrom}'")

            p_start, p_end, q_start, q_end = FacetsMeta.chr_arms[key]

            # interval entirely in p
            if end <= p_end:
                arm = "p"
            # entirely in q
            elif start >= q_start:
                arm = "q"
            # spans centromere
            else:
                arm = "pq"

            return f"{chrom}{arm}"

    def setVerbose(self, doVerbose):
        """Print per-sample diagnostics (skipped fits, unreadable manifests) while building."""
        self.run_verbose = bool(doVerbose)

    def _vprint(self, message, color=bcolors.WARNING):
        if self.run_verbose:
            print(color + message + bcolors.ENDC)

    #This function will build the FacetsMeta object. It handles executing parseClinicalSample() and data persistence.
    def buildFacetsMeta(self):
        #Making a meta is going to be the first thing in FACETS API, so this is where we can print the logo.
        self.printLogo()

        #Fill in most of these structures with the parseClinicalSample function.
        #If we are using data persistance, we want to load that from an existing structured file or save the initial parse to an existing file.
        #That way we don't need to process the clinical sample file each run.
        if self.persist_data != "no":
            if os.path.isfile(self.persist_data):
                print("\tLoading previously processed FacetsMeta object: " + self.persist_data)
                with open(self.persist_data, 'rb') as storedMeta:
                    savedInstance = pickle.load(storedMeta)
                saved = dict(getattr(savedInstance, "__dict__", {}))
                version = saved.get("_persist_version")
                if version != FacetsMeta.PERSIST_VERSION:
                    print(bcolors.WARNING + "\t\tWarning: this persisted FacetsMeta was written by an older facetsAPI "
                          "(format %s, current %s). Fit-selection settings are taken from this object, not from the file; "
                          "rebuild the persisted file to refresh it." % (version, FacetsMeta.PERSIST_VERSION) + bcolors.ENDC)
                    for k in ("selectSingleRun", "use_unreviewed_defaults", "fit_policy", "_fit_policy_explicit",
                              "fit_class", "hisens_vs_purity", "_persist_version"):
                        saved.pop(k, None)
                for k in ("persist_data", "run_verbose"):
                    saved.pop(k, None)
                for k, v in saved.items():
                    setattr(self, k, v)
            else:
                self.parseClinicalSample()
                print("\tStoring this FacetsMeta object for future loading: " + self.persist_data)
                with open(self.persist_data, 'wb') as dataToStore:
                    pickle.dump(self, dataToStore)
        else:
            self.parseClinicalSample()

    #This function will set this FacetsMeta object to select a single run per sample. If allowDefaults is set to True,
    # default fits will be selected when reviewed/autoQC fits are not available.  If false, only reviewed/autoQC fits will be selected.
    #Use setFitPolicy() for finer control (reviewed only, autoqc only, default only, tier order).
    def setSingleRunPerSample(self, doSetSingle, allowDefaults=False):
        if not isinstance(doSetSingle, bool):
            raise FacetsConfigError("FacetsMeta.setSingleRunPerSample(): doSetSingle expected True or False.")
        if not isinstance(allowDefaults, bool):
            raise FacetsConfigError("FacetsMeta.setSingleRunPerSample(): allowDefaults expected True or False.")
        if doSetSingle is False and allowDefaults is True:
            raise FacetsConfigError("FacetsMeta.setSingleRunPerSample(): allowDefaults can only be True when doSetSingle is also True.")

        self.selectSingleRun         = doSetSingle
        self.use_unreviewed_defaults = allowDefaults
        if doSetSingle:
            if self._fit_policy_explicit:
                #Keep the caller's policy, only toggle the default tier.
                self.fit_policy = self.fit_policy.with_defaults(allowDefaults)
            elif allowDefaults:
                self.fit_policy = FitPolicy((TIER_REVIEWED, TIER_AUTOQC, TIER_DEFAULT))
            else:
                self.fit_policy = FitPolicy((TIER_REVIEWED, TIER_AUTOQC))

    #Set the fit-selection policy used in single-run mode. Examples:
    #   FitPolicy()                       reviewed > autoqc > default (the default)
    #   FitPolicy.reviewed_only()         only human-reviewed best/acceptable fits
    #   FitPolicy.autoqc_only()           only autoQC best/pass fits
    #   FitPolicy.default_only()          only the facets default fit
    #   FitPolicy(("autoqc", "reviewed")) autoQC preferred over human review
    def setFitPolicy(self, policy):
        if not isinstance(policy, FitPolicy):
            raise FacetsConfigError("FacetsMeta.setFitPolicy(): policy must be a FitPolicy.")
        self.fit_policy              = policy
        self._fit_policy_explicit    = True
        self.use_unreviewed_defaults = policy.allows_defaults

    #This function accepts a file with a single sample ID per line and populates this objects
    #samples_from_file list.  If this list is populated and build_from_file_listing is true,
    #When parseClinicalData runs, it will only include samples listed in the provided file.
    #Entries may be tumor sample ids (P-0000001-T01-IM6) or full pair ids; blank and '#' lines are ignored.
    def selectSamplesFromFile(self, infile):
        try:
            with open(infile) as fp:
                for line in fp:
                    entry = line.strip()
                    if entry and not entry.startswith("#"):
                        self.samples_from_file.append(entry)
        except OSError as e:
            raise FacetsRepositoryError("FacetsMeta.selectSamplesFromFile(): cannot read %s: %s" % (infile, e)) from e
        self.build_from_file_listing = True

    #Reduce any master_file_dict key / run id to the 17-character tumor sample id.
    @staticmethod
    def shortId(key):
        return short_sample_id(key)

    #The FacetsSample id (<tumor>_<normal> pair id) that a master_file_dict key belongs to.
    def pairIdForKey(self, key):
        base = self.long_id_map.get(short_sample_id(key))
        if base:
            return pair_id(base[0])
        return pair_id(key)

    #Simple method for printing a FACETS API logo.
    @staticmethod
    def printLogo():
        print("~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~") 
        print("|"+bcolors.OKBLUE+r" ______   ______     ______     ______     ______   ______    "+bcolors.OKCYAN+r"    ______     ______   __    "+bcolors.ENDC+"|") 
        print("|"+bcolors.OKBLUE+r"/\  ___\ /\  __ \   /\  ___\   /\  ___\   /\__  _\ /\  ___\   "+bcolors.OKCYAN+r"   /\  __ \   /\  == \ /\ \   "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+r"\ \  __\ \ \  __ \  \ \ \____  \ \  __\   \/_/\ \/ \ \___  \  "+bcolors.OKCYAN+r"   \ \  __ \  \ \  _-/ \ \ \  "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+r" \ \_\    \ \_\ \_\  \ \_____\  \ \_____\    \ \_\  \/\_____\ "+bcolors.OKCYAN+r"    \ \_\ \_\  \ \_\    \ \_\ "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+r"  \/_/     \/_/\/_/   \/_____/   \/_____/     \/_/   \/_____/  "+bcolors.OKCYAN+r"    \/_/\/_/   \/_/     \/_/ "+bcolors.ENDC+"|")
        print("~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~")
        
    ######################
    # parseClinicalSample:  This function will accept a clinical sample file and
    #                         scan each sample's corresponding facets directory.
    #                         Each sample directory is inspected with layout.discover_sample(): a standard
    #                         sample has facets_review.manifest at its root; a 2n pair has clinical/ and
    #                         research/ subdirectories that each carry their own manifest. The manifest(s)
    #                         are read to identify the fit(s) to use, and a dictionary mapping run keys to
    #                         their files is built in master_file_dict. Fit directories are resolved from
    #                         the DISCOVERED directory, never from the manifest's absolute 'path' column.
    #                         This function also builds dictionaries mapping sample id to cancer type,
    #                         patient id, etc. from the clinical sample file.
    ######################
    COUNTER_KEYS = ("samples_considered", "not_found", "multi_normal", "missing_manifest", "missing_class",
                    "bad_manifest", "no_fit_selected", "missing_files", "skipped_ultra", "skipped_not_selected", "duplicate_rows",
                    "standard_samples", "twon_samples", "runs_registered",
                    "Best", "Acceptable", "AutoQC_Best", "AutoQC_Pass", "Default", "Unreviewed")
    CLINICAL_REQUIRED_COLUMNS = ("SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "TUMOR_PURITY",
                                 "ONCOTREE_CODE", "CVR_TMB_SCORE", "MSI_SCORE")

    def _sampleFilterSet(self):
        return {str(s).strip()[0:17] for s in self.samples_from_file if str(s).strip()}

    def _readClinicalFile(self):
        try:
            clinical_df = pd.read_csv(self.clinical_sample_file, sep="\t", low_memory=False)
        except OSError as e:
            raise FacetsRepositoryError("Cannot read clinical sample file %s: %s" % (self.clinical_sample_file, e)) from e
        except (pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            raise FacetsRepositoryError("Malformed clinical sample file %s: %s" % (self.clinical_sample_file, e)) from e
        missing = [c for c in FacetsMeta.CLINICAL_REQUIRED_COLUMNS if c not in clinical_df.columns]
        if missing:
            raise FacetsRepositoryError("Clinical sample file %s lacks column(s): %s"
                                        % (self.clinical_sample_file, ", ".join(missing)))
        return clinical_df

    def _registerRun(self, key, fc, cand, files, counts):
        """Store one fit under ``key`` in master_file_dict (MetaDictMap layout) and count it."""
        selected_fit_dir = _fit_dir(fc.class_dir, cand.fit_name) + "/"
        self.master_file_dict[key] = [files["out"], files["cncf"], files["qc"], fc.facets_qc_path, selected_fit_dir,
                                      files["gene"], files["adjseg"], files["ccf"], files["nonsignedout"],
                                      files["unadjseg"], fc.manifest_path, cand.fit_status, fc.class_dir,
                                      files["rdata"], fc.name, cand.tier]
        self.fit_map[key] = [cand.label, selected_fit_dir]
        counts[cand.label] = counts.get(cand.label, 0) + 1
        counts["runs_registered"] += 1
        return selected_fit_dir

    def _selectSingleFit(self, id, id_with_normal, fc, manifest, multi_class, counts):
        cand, files, tried = self.fit_policy.select_fit(manifest, self.hisens_vs_purity, fc.class_dir)
        if cand is None:
            if tried:
                counts["missing_files"] += 1
                self._vprint("\t\tSkipping %s (%s): candidate fit(s) %s lack %s files."
                             % (id_with_normal, fc.name, [t.fit_name for t in tried], self.hisens_vs_purity))
            else:
                counts["no_fit_selected"] += 1
                self._vprint("\t\tSkipping %s (%s): no fit matches policy tiers %s."
                             % (id_with_normal, fc.name, list(self.fit_policy.tiers)))
            return
        key = id if not multi_class else id + "#" + fc.name
        self.long_id_map.setdefault(id, [id_with_normal])
        self._registerRun(key, fc, cand, files, counts)
        if key != id and id not in self.fit_map:
            self.fit_map[id] = list(self.fit_map[key])

    def _collectAllFits(self, id, id_with_normal, fc, manifest, counts):
        mode = self.hisens_vs_purity
        #Label the policy-preferred fit for fit_map[tumor id] (what writeReport's fitType column shows).
        pref, pref_files, _tried = self.fit_policy.select_fit(manifest, mode, fc.class_dir)
        if pref is not None and id not in self.fit_map:
            self.fit_map[id] = [pref.label, _fit_dir(fc.class_dir, pref.fit_name) + "/"]
        rows, n_dup = self.fit_policy.representative_rows(manifest)
        counts["duplicate_rows"] += n_dup
        for idx, row in rows.iterrows():
            fit_name = str(row["fit_name"])
            status   = str(row["review_status"])
            files, missing = locate_fit_files(_fit_dir(fc.class_dir, fit_name), mode)
            if files is None:
                counts["missing_files"] += 1
                self._vprint("\t\tSkipping fit %s of %s (%s): missing %s." % (fit_name, id_with_normal, fc.name, missing))
                continue
            cand = FitCandidate(fit_name, status, str(row["reviewed_by"]), tier_of(status, fit_name),
                                label_of(status, fit_name), int(idx))
            fit_label = fit_name if not fc.is_2n else fc.name + "/" + fit_name
            long_id = id_with_normal + "#" + fit_label
            self.long_id_map.setdefault(id, []).append(long_id)
            self._registerRun(long_id, fc, cand, files, counts)

    def parseClinicalSample(self):
        counts = {k: 0 for k in FacetsMeta.COUNTER_KEYS}
        self.build_counts = counts

        print("\tProcessing clinical data...")
        if self.facets_repo_path == "":
            raise FacetsRepositoryError("FacetsMeta.parseClinicalSample(): Facets repo path is not supplied.")
        repo = self.facets_repo_path
        if not os.path.isdir(repo):
            raise FacetsRepositoryError("FacetsMeta.parseClinicalSample(): Facets repo path is not a directory: %s" % repo)

        #If there is no clinical sample file,
        #we just want to build the FacetsMeta object by browsing the directory stucture.
        if self.clinical_sample_file == "":
            print(bcolors.WARNING + "\t\tWarning: No data clinical sample file was provided. ")
            print("\t\t\tUsing FACETS directory structure to build FacetsMeta object: " + repo + " " + bcolors.ENDC)
            target_ids = []
            for item in sorted(os.listdir(repo)):
                sub = os.path.join(repo, item)
                if os.path.isdir(sub):
                    for samplefolder in sorted(os.listdir(sub)):
                        target_ids.append(samplefolder)
        else:
            clinical_df = self._readClinicalFile()

            #Extract DMP id's for everything in our clinical sample file.
            print("\t\tIdentifying samples from impact FACETS repository.")

            target_ids           = clinical_df['SAMPLE_ID'].astype(str).tolist()
            patient_ids          = clinical_df['PATIENT_ID'].tolist()
            cancer_types         = clinical_df['CANCER_TYPE'].tolist()
            cancer_detail_types  = clinical_df['CANCER_TYPE_DETAILED'].tolist()
            clin_purities        = clinical_df['TUMOR_PURITY'].tolist()
            oncotree_codes       = clinical_df['ONCOTREE_CODE'].tolist()
            tmb_scores           = clinical_df['CVR_TMB_SCORE'].tolist()
            msi_scores           = clinical_df['MSI_SCORE'].tolist()

            #Build sample id maps associated with cancer type and patient ID.
            print("\t\tCreating data to id maps.")
            for i in range(len(target_ids)):
                self.patient_id_map[target_ids[i]]         = patient_ids[i]
                self.cancer_type_map[target_ids[i]]        = cancer_types[i]
                self.cancer_type_detail_map[target_ids[i]] = cancer_detail_types[i]
                self.clinical_purity_map[target_ids[i]]    = clin_purities[i]
                self.onkotree_code_map[target_ids[i]]      = oncotree_codes[i]
                self.cvr_tmb_score_map[target_ids[i]]      = tmb_scores[i]
                self.msi_score_map[target_ids[i]]          = msi_scores[i]

        #For each DMP id, confirm existence and build a set of relevant directories.
        if self.selectSingleRun:
            print("\t\tBuilding directory map targeting one run per sample (policy tiers: %s)." % " > ".join(self.fit_policy.tiers))
        else:
            print("\t\tBuilding directory map including all runs for each sample.")
        print("\t\tReading %s runs; 2n fit class: %s." % (self.hisens_vs_purity, self.fit_class))

        sample_filter = None
        if self.build_from_file_listing:
            sample_filter = self._sampleFilterSet()
            print("\t\tSelecting samples only from the provided input sample list (%d ids)." % len(sample_filter))

        for id in target_ids:
            id = str(id)
            cur_tumor_id = id[0:17]

            #If we are parsing data from a file, only include the samples we've already read in.
            if sample_filter is not None and cur_tumor_id not in sample_filter:
                continue
            counts["samples_considered"] += 1

            cur_id_dirs = sorted(d for d in glob.glob(os.path.join(repo, id[0:7], id + "*")) if os.path.isdir(d))
            if not cur_id_dirs:
                counts["not_found"] += 1
                continue
            #Keep track of how many samples had the same tumor/id, but multiple normals. The last one (sorted) is used.
            if len(cur_id_dirs) > 1:
                counts["multi_normal"] += 1
            sample_dir     = cur_id_dirs[-1].rstrip("/")
            id_with_normal = os.path.basename(sample_dir)

            classes = discover_sample(sample_dir)
            if not classes:
                counts["missing_manifest"] += 1
                continue
            if classes[0].is_2n:
                counts["twon_samples"] += 1
            else:
                counts["standard_samples"] += 1
            classes = select_classes(classes, self.fit_class)
            if not classes:
                counts["missing_class"] += 1
                self._vprint("\t\tSkipping %s: 2n pair has no %s class." % (id_with_normal, self.fit_class))
                continue
            multi_class = len(classes) > 1

            for fc in classes:
                try:
                    manifest = read_manifest(fc.manifest_path)
                except ManifestError as e:
                    counts["bad_manifest"] += 1
                    self.failed_samples.append(id_with_normal)
                    self._vprint("\t\t" + str(e), bcolors.FAIL)
                    continue
                counts["skipped_not_selected"] += int((manifest["fit_name"] == NOT_SELECTED).sum())
                if not self.fit_policy.include_ultra_fits:
                    counts["skipped_ultra"] += int(manifest["fit_name"].map(is_ultra).sum())

                if self.selectSingleRun:
                    self._selectSingleFit(id, id_with_normal, fc, manifest, multi_class, counts)
                else:
                    self._collectAllFits(id, id_with_normal, fc, manifest, counts)

        if self.selectSingleRun:
            print("\tSelecting single FacetsRun per FacetsSample.")
            print("\t\tIdentified " + str(counts["Best"]) + " reviewed best fits.")
            print("\t\tIdentified " + str(counts["Acceptable"]) + " reviewed acceptable fits.")
            print("\t\tIdentified " + str(counts["AutoQC_Best"]) + " autoQC best fits.")
            print("\t\tIdentified " + str(counts["AutoQC_Pass"]) + " autoQC pass fits.")
            print("\t\tIdentified " + str(counts["Default"]) + " default fits.")
            print("\t\t" + str(counts["no_fit_selected"]) + " samples had no fit matching the policy.")
        else:
            print("\tSelecting all FacetsRuns for each FacetsSample.")
        print("\tLayouts: %d standard samples, %d 2n pairs (fit_class=%s)."
              % (counts["standard_samples"], counts["twon_samples"], self.fit_class))
        print("\tIdentified " + str(counts["multi_normal"]) + " instances of tumor having multiple normals.")
        print("\tSkipped %d ids with no directory, %d samples with no manifest, %d 2n pairs missing the requested class, %d unreadable manifests."
              % (counts["not_found"], counts["missing_manifest"], counts["missing_class"], counts["bad_manifest"]))
        print("\tSkipped " + str(counts["missing_files"]) + " fits due to missing files.")
        if counts["skipped_ultra"]:
            print("\tExcluded " + str(counts["skipped_ultra"]) + " ultra fits (FitPolicy.include_ultra_fits is False).")
        if counts["duplicate_rows"]:
            print("\tCollapsed " + str(counts["duplicate_rows"]) + " duplicate manifest rows (same fit listed more than once; highest-ranked row kept).")
        print("\t\tTotal runs extracted: " + str(len(self.master_file_dict)))


######################
# FacetsDataset:    This class represents a set of FacetsRuns.  It includes functions for selecting, filtering,
#                   writing output, and modifying/creating file structures.
######################
class FacetsDataset:
    def __init__(self, facets_meta):
        self.sampleList = {}
        self.runList = []
        self.numSamples = 0
        self.ref_facetsMeta = facets_meta

        self.filterCancerType = False
        self.selectedCancerTypes = []
        self.filterCancerTypeDetail = False
        self.selectedCancerDetailTypes = []
        self.filterPurity = False
        self.minPurity = -1.0
        self.maxPurity = -1.0
        self.filterClinicalPurity = False
        self.minClinPurity = -1.0
        self.maxClinPurity = -1.0
        self.filterOnkoCode = False
        self.selectedOnkoCodes = []
        self.filterPloidy = False
        self.minPloidy    = -1.0
        self.maxPloidy    = -1.0
        self.filterDipLogR = False
        self.minDipLogR    = -1.0
        self.maxDipLogR    = -1.0
        self.filterCval = False
        self.minC       = -1.0
        self.maxC       = -1.0
        self.filterWGD = False
        self.keepWGD   = False
        self.filterFGA = False
        self.minFGA    = -1.0
        self.maxFGA    = -1.0
        self.filterFracLoH = False
        self.minLoh        = -1.0
        self.maxLoh        = -1.0
        self.filterFacetsQC = False
        self.keepQCPass = False
        self.filterTMB = False
        self.minTMB    = -1.0
        self.maxTMB    = -1.0
        self.filterMSI = False
        self.minMSI    = -1.0
        self.maxMSI    = -1.0
        self.num_build_failed = 0
        self.num_file_failed = 0

    #This function will write a report for the samples in the current dataset to the given outfile.
    def writeReport(self, outfile_name):
        print("\t\t\tWriting report to: " + outfile_name)
        with open(outfile_name, 'w') as outfile:
            outfile.write("ID\tfitDir\tfitType\tCancer Type\tCancer Type Detail\tFacets Purity\tClinical Purity\tOnkoCode\tPloidy\tdipLogR\tcVal\tWGD\tFGA\tfrac_loh\tfacets_qc_pass\ttmb\tmsi" + "\n")
            found_list = []
            for cur_run in self.runList:
                #fit_map is keyed by the run key (single-run) or long id (all fits); fall back to the tumor id.
                fit_map = self.ref_facetsMeta.fit_map
                cur_fit_item = fit_map.get(cur_run.id) or fit_map.get(short_sample_id(cur_run.id))

                outfile.write(str(cur_run.id) + "\t")
                outfile.write(str(cur_run.fitDir) + "\t")
                if cur_fit_item is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_fit_item[0]) + "\t")
                if cur_run.cancerType is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.cancerType) + "\t")
                if cur_run.cancerTypeDetail is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.cancerTypeDetail) + "\t")
                if cur_run.purity is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.purity) + "\t")
                if cur_run.clinicalPurity is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.clinicalPurity) + "\t")
                if cur_run.onkoCode is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.onkoCode) + "\t")
                outfile.write(str(cur_run.ploidy) + "\t")
                outfile.write(str(cur_run.dipLogR) + "\t")
                outfile.write(str(cur_run.cval) + "\t")
                outfile.write(str(cur_run.wgd) + "\t")
                outfile.write(str(cur_run.fga) + "\t")
                outfile.write(str(cur_run.frac_loh) + "\t")
                outfile.write(str(cur_run.facets_qc) + "\t")
                if cur_run.tmb is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.tmb) + "\t")
                if cur_run.msi is None:
                    outfile.write("NA" + "\t")
                else:
                    outfile.write(str(cur_run.msi) + "\t")
                outfile.write("\n")
    #This function will make backup copies (<file>.<date>.bak) of the ccf.maf file of every run in the provided FacetsMeta object.
    def backupMafs(self, ref_meta):
        if not isinstance(ref_meta, FacetsMeta):
            raise FacetsConfigError("FacetsDataset.backupMafs(): ref_meta must be of type FacetsMeta.")
        print("Backing up ccf.maf files.")
        done = set()
        for item in ref_meta.master_file_dict:
            cur_maf = ref_meta.master_file_dict.get(item)[MetaDictMap.CCF_MAF_FILE]
            if not cur_maf or cur_maf in done or not os.path.isfile(cur_maf):
                continue
            done.add(cur_maf)
            shutil.copy2(cur_maf, cur_maf + "." + str(date.today()) + ".bak")
        return len(done)

    #This function outputs a text file containing information that summarizes a facets dataset
    #This includes metrics like number of best/acceptable fits, number of failed samples that couldn't build, number of cancer types, etc
    #The output contains the date and time.  
    def writeDatasetSummary(self, outdir="."):
        status_dict         = {}
        cancer_type_dict    = {}
        oncocode_dict       = {}
        msi_count           = 0
        msi_li              = []
        wgd_count           = 0
        purity_li           = []
        bestfit_count       = 0
        acceptablefit_count = 0
        default_count       = 0

        for run in self.runList:
            if run.review_status not in status_dict:
                status_dict[run.review_status] = 0
            status_dict[run.review_status] += 1

            if run.fit_name not in status_dict:
                status_dict[run.fit_name] = 0
            status_dict[run.fit_name] += 1
            
            if run.fit_name == "default":
                default_count +=1
            if "best" in run.review_status: 
                bestfit_count +=1
            if "acceptable" in run.review_status: 
                acceptablefit_count +=1

            if run.cancerType is not None:
                curCancerType = run.cancerType
                if curCancerType not in cancer_type_dict:
                    cancer_type_dict[curCancerType] = 0
                cancer_type_dict[curCancerType] += 1
            
            if run.onkoCode is not None:
                curOnkoCode = run.onkoCode
                if curOnkoCode not in oncocode_dict:
                    oncocode_dict[curOnkoCode] = 0
                oncocode_dict[curOnkoCode] += 1

            if run.msi is not None:
                cur_msi = run.msi
                if cur_msi>=1.5:
                    msi_count +=1
                msi_li.append(cur_msi)

            if run.wgd:
                wgd_count+=1
            
            purity_li.append(run.purity)
        
        msi_median = statistics.median(msi_li) if msi_li else "NA"
        purity_vals = [v for v in (_to_float(p) for p in purity_li) if v is not None]
        purity_median = statistics.median(purity_vals) if purity_vals else "NA"

        date = datetime.now()
        datestr = str(date.isoformat(sep = "_",timespec='minutes')).replace(":","-")
        outstr = ""
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "FacetsDataset_summary_{datest}.tsv".format(datest= datestr)),'w') as outf:
            outstr+= "FacetsDataset_summary_{datest}".format(datest= datestr) + '\n'
            outstr+= "Number of Samples" + '\t' + str(self.numSamples) + '\n'
            outstr+= "Number of files failed" + '\t' + str(self.num_file_failed )+ '\n'
            outstr+= "Number of failed sample build" + '\t' + str(self.num_build_failed) + '\n'

            for key in status_dict:
                outstr+= key+" count\t"+str(status_dict[key])+'\n'

            outstr+= "Cancer Type Number"+"\t"+str(len(cancer_type_dict))+'\n'
                
            for key in cancer_type_dict:
                outstr+= key+" count\t"+str(cancer_type_dict[key])+'\n'

            outstr+= "Onco Code Number"+"\t"+str(len(oncocode_dict))+'\n'
                
            for key in oncocode_dict:
                outstr+= key+" count\t"+str(oncocode_dict[key])+'\n'

            outstr += "Best fit Count"+"\t"+str(bestfit_count)+'\n'
            outstr += "Acceptable fit Count"+"\t"+str(acceptablefit_count)+'\n'
            outstr += "Default fit Count"+"\t"+str(default_count)+'\n'
            outstr += "MSI score median "+"\t"+str(msi_median)+'\n'
            outstr += "MSI above 1.5"+"\t"+str(msi_count)+'\n'
            outstr += "Purity median"+"\t"+str(purity_median)+'\n'
            outstr += "WGD count"+"\t"+str(wgd_count)+'\n'
            
            outf.write(outstr)


    ###################### 
    # createHistogram:  This function outputs a histogram png using the variable requested from facetsrun
    #                   It defaults to 10 automated bins for numerical variables but that can be changed. 
    #                   For strings like cancerType or boolean values like WGD it uses however many bins are needed.
    #                   Available variables are as follows:   'purity', 'clinicalPurity',  'ploidy', 'dipLogR', 'cval', 'fga', 'frac_loh', 'tmb' , 'msi'
    #                   'cancerType', 'oncoCode','wgd', 'review_status'
    ######################
    NUMERIC_HISTOGRAM_VARIABLES     = ('purity', 'clinicalPurity', 'ploidy', 'dipLogR', 'cval', 'fga', 'frac_loh', 'tmb', 'msi')
    CATEGORICAL_HISTOGRAM_VARIABLES = ('cancerType', 'onkoCode', 'oncoCode', 'wgd', 'review_status', 'fit_class', 'fit_name')

    def createHistogram(self, variable, inbins=10, outdir="."):
        """Write FacetsDataset_histogram_<variable>_<datetime>.png into ``outdir`` and return its path."""
        attr = 'onkoCode' if variable == 'oncoCode' else variable   # accept the documented spelling too
        date = datetime.now()
        datestr = str(date.isoformat(sep="_", timespec='minutes')).replace(":", "-")
        os.makedirs(outdir, exist_ok=True)
        out_path = os.path.join(outdir, 'FacetsDataset_histogram_{var}_{datetim}.png'.format(var=variable, datetim=datestr))
        plt = _pyplot()

        if variable in FacetsDataset.NUMERIC_HISTOGRAM_VARIABLES:
            variableli = [v for v in (_to_float(getattr(run, attr, None)) for run in self.runList) if v is not None]
            n, bins, patches = plt.hist(x=variableli, bins=inbins, color='#0504aa',
                                        alpha=0.7, rwidth=0.85, orientation="horizontal")
            plt.ylabel('Value')
            plt.xlabel('Frequency')
            plt.title("Histogram of {var} in constructed FacetsDataset ({datetim})".format(var=variable, datetim=datestr))
            plt.yticks(bins)
        elif variable in FacetsDataset.CATEGORICAL_HISTOGRAM_VARIABLES:
            count_dict = {}
            for run in self.runList:
                label = str(getattr(run, attr, None))
                count_dict[label] = count_dict.get(label, 0) + 1
            plt.barh(list(count_dict.keys()), list(count_dict.values()), color='g')
            plt.xlabel('Value')
            plt.ylabel('Frequency')
            plt.title("Histogram of {var} in constructed FacetsDataset ({datetim})".format(var=variable, datetim=datestr))
            plt.autoscale()
        else:
            raise FacetsConfigError("FacetsDataset.createHistogram(): unknown variable %r; choose from %s"
                                    % (variable, FacetsDataset.NUMERIC_HISTOGRAM_VARIABLES + FacetsDataset.CATEGORICAL_HISTOGRAM_VARIABLES))
        plt.savefig(out_path, bbox_inches="tight")
        plt.clf()
        return out_path


    ######################
    # copyDatasetToFolder:  This function will create a copy of a the folders in the facets dataset object
    #                       to a specified output location. Directory structure will be preserved, /all/P-12345/P-12345_T01...
    #                       Note that this function will not overwrite existing files.
    ######################     
    def copyDatasetToFolder(self, output_folder):
        """Copy every sample directory in this dataset to <output_folder>/all/<P-00xxx>/<pair>/ (no overwrite).

        For 2n pairs the whole pair directory (clinical/ and research/) is copied.
        """
        if os.path.abspath(output_folder) == os.path.abspath(self.ref_facetsMeta.facets_repo_path):
            raise FacetsConfigError("FacetsDataset.copyDatasetToFolder(): destination path cannot be the same as the source path.")

        target_base_dir = os.path.join(output_folder, "all")
        for d in (output_folder, target_base_dir):
            if not os.path.exists(d):
                os.makedirs(d, exist_ok=True)
                os.chmod(d, 0o775)

        copied = 0
        for key in self.sampleList:
            cur_run = self.sampleList.get(key).runs[0]  # all runs of a sample share the sample directory
            base_dir = self.ref_facetsMeta.master_file_dict.get(cur_run.id, [None] * 13)[MetaDictMap.SAMPLE_BASE_DIR]
            if not base_dir:
                base_dir = os.path.dirname(os.path.dirname(cur_run.fitDir.rstrip("/")))
            base_dir = base_dir.rstrip("/")
            if cur_run.fit_class != "standard":            # 2n: step up from the class dir to the pair dir
                base_dir = os.path.dirname(base_dir)
            long_dir_str  = os.path.basename(base_dir)                   # <tumor>_<normal>
            short_dir_str = os.path.basename(os.path.dirname(base_dir))  # P-00xxx
            new_short_parent_dir = os.path.join(target_base_dir, short_dir_str)
            if not os.path.exists(new_short_parent_dir):
                os.makedirs(new_short_parent_dir, exist_ok=True)
                os.chmod(new_short_parent_dir, 0o775)
            target_copy_path = os.path.join(new_short_parent_dir, long_dir_str)
            if not os.path.exists(target_copy_path):
                shutil.copytree(base_dir, target_copy_path, symlinks=True)
                copied += 1
        print("\t\tCopied %d sample directories to %s" % (copied, target_base_dir))
        return copied


    #This function will calculate and write calculated purity data for each sample.  If use_base_cf is true, cf will be used
    #rather than cf.em.  
    def writePurityCFs(self, outfile_path, use_base_cf=False, use_long_id=True):
        with open(outfile_path, 'w') as outfile:
            headerString = "ID\tPurity\n"
            outfile.write(headerString)
            for curSample in self.runList:
                curPurity = curSample.calculatePurityByCF(use_base_cf)
                if use_long_id:
                    outfile.write(self.ref_facetsMeta.pairIdForKey(curSample.id) + "\t" + str(curPurity) + "\n")
                else:
                    outfile.write(curSample.id + "\t" + str(curPurity) + "\n")

    #This function will run alteration analysis on the FacetsDataset.
    def runAlterationAnalysis(self):
        for sample in self.sampleList:
            cur_sample = self.sampleList.get(sample)
            for cur_run in cur_sample.runs:
                cur_run.defineArms()
                cur_run.defineAllLohAndGains()
                cur_run.defineCFLevels()
                cur_run.defineSampleLevelAlterationCall(FacetsRun.alteration_cov_min, FacetsRun.alteration_arm_num)
    #This function will write an output file containing arm level data and alteration calls.
    #Output file will be id, cancer_type, cancer_type_detail, purity, onkoCode, tmb, msi, wgd, fga, frac_loh, isHypoploidy, isGainHeavy
    #                        1p_alteration, 1q_alteration...22q_alteration, 
    #                        1p_percentAltered, 1q_percentAltered...22q_percentAltered
    #                        1p_cf_value, 1q_cf_value...22q_cf_value, 
    #                        1p_cf_armCov, 1q_cf_armCov...22q_cf_value,
    #                        1p_level, 1q_level...22q_level
    def writeAlterationData(self, outfile_path):
        with open(outfile_path, 'w') as outfile:
            headerString = "ID\t"
            headerString += "CancerType\t"
            headerString += "CancerTypeDetail\t"
            headerString += "Purity\t"
            headerString += "OnkoCode\t"
            headerString += "TMB\t"
            headerString += "MSI\t"
            headerString += "WGD\t"
            headerString += "FGA\t"
            headerString += "FracLoH\t"
            headerString += "isHypoploidy\t"
            headerString += "isGainHeavy\t"

            headerString += "1p_alteration\t"
            headerString += "1q_alteration\t"
            headerString += "2p_alteration\t"
            headerString += "2q_alteration\t"
            headerString += "3p_alteration\t"
            headerString += "3q_alteration\t"
            headerString += "4p_alteration\t"
            headerString += "4q_alteration\t"
            headerString += "5p_alteration\t"
            headerString += "5q_alteration\t"
            headerString += "6p_alteration\t"
            headerString += "6q_alteration\t"
            headerString += "7p_alteration\t"
            headerString += "7q_alteration\t"
            headerString += "8p_alteration\t"
            headerString += "8q_alteration\t"
            headerString += "9p_alteration\t"
            headerString += "9q_alteration\t"
            headerString += "10p_alteration\t"
            headerString += "10q_alteration\t"
            headerString += "11p_alteration\t"
            headerString += "11q_alteration\t"
            headerString += "12p_alteration\t"
            headerString += "12q_alteration\t"
            headerString += "13p_alteration\t"
            headerString += "13q_alteration\t"
            headerString += "14p_alteration\t"
            headerString += "14q_alteration\t"
            headerString += "15p_alteration\t"
            headerString += "15q_alteration\t"
            headerString += "16p_alteration\t"
            headerString += "16q_alteration\t"
            headerString += "17p_alteration\t"
            headerString += "17q_alteration\t"
            headerString += "18p_alteration\t"
            headerString += "18q_alteration\t"
            headerString += "19p_alteration\t"
            headerString += "19q_alteration\t"
            headerString += "20p_alteration\t"
            headerString += "20q_alteration\t"
            headerString += "21p_alteration\t"
            headerString += "21q_alteration\t"
            headerString += "22p_alteration\t"
            headerString += "22q_alteration\t"

            headerString += "1p_alteredPercent\t"
            headerString += "1q_alteredPercent\t"
            headerString += "2p_alteredPercent\t"
            headerString += "2q_alteredPercent\t"
            headerString += "3p_alteredPercent\t"
            headerString += "3q_alteredPercent\t"
            headerString += "4p_alteredPercent\t"
            headerString += "4q_alteredPercent\t"
            headerString += "5p_alteredPercent\t"
            headerString += "5q_alteredPercent\t"
            headerString += "6p_alteredPercent\t"
            headerString += "6q_alteredPercent\t"
            headerString += "7p_alteredPercent\t"
            headerString += "7q_alteredPercent\t"
            headerString += "8p_alteredPercent\t"
            headerString += "8q_alteredPercent\t"
            headerString += "9p_alteredPercent\t"
            headerString += "9q_alteredPercent\t"
            headerString += "10p_alteredPercent\t"
            headerString += "10q_alteredPercent\t"
            headerString += "11p_alteredPercent\t"
            headerString += "11q_alteredPercent\t"
            headerString += "12p_alteredPercent\t"
            headerString += "12q_alteredPercent\t"
            headerString += "13p_alteredPercent\t"
            headerString += "13q_alteredPercent\t"
            headerString += "14p_alteredPercent\t"
            headerString += "14q_alteredPercent\t"
            headerString += "15p_alteredPercent\t"
            headerString += "15q_alteredPercent\t"
            headerString += "16p_alteredPercent\t"
            headerString += "16q_alteredPercent\t"
            headerString += "17p_alteredPercent\t"
            headerString += "17q_alteredPercent\t"
            headerString += "18p_alteredPercent\t"
            headerString += "18q_alteredPercent\t"
            headerString += "19p_alteredPercent\t"
            headerString += "19q_alteredPercent\t"
            headerString += "20p_alteredPercent\t"
            headerString += "20q_alteredPercent\t"
            headerString += "21p_alteredPercent\t"
            headerString += "21q_alteredPercent\t"
            headerString += "22p_alteredPercent\t"
            headerString += "22q_alteredPercent\t"

            headerString += "1p_cf\t"
            headerString += "1q_cf\t"
            headerString += "2p_cf\t"
            headerString += "2q_cf\t"
            headerString += "3p_cf\t"
            headerString += "3q_cf\t"
            headerString += "4p_cf\t"
            headerString += "4q_cf\t"
            headerString += "5p_cf\t"
            headerString += "5q_cf\t"
            headerString += "6p_cf\t"
            headerString += "6q_cf\t"
            headerString += "7p_cf\t"
            headerString += "7q_cf\t"
            headerString += "8p_cf\t"
            headerString += "8q_cf\t"
            headerString += "9p_cf\t"
            headerString += "9q_cf\t"
            headerString += "10p_cf\t"
            headerString += "10q_cf\t"
            headerString += "11p_cf\t"
            headerString += "11q_cf\t"
            headerString += "12p_cf\t"
            headerString += "12q_cf\t"
            headerString += "13p_cf\t"
            headerString += "13q_cf\t"
            headerString += "14p_cf\t"
            headerString += "14q_cf\t"
            headerString += "15p_cf\t"
            headerString += "15q_cf\t"
            headerString += "16p_cf\t"
            headerString += "16q_cf\t"
            headerString += "17p_cf\t"
            headerString += "17q_cf\t"
            headerString += "18p_cf\t"
            headerString += "18q_cf\t"
            headerString += "19p_cf\t"
            headerString += "19q_cf\t"
            headerString += "20p_cf\t"
            headerString += "20q_cf\t"
            headerString += "21p_cf\t"
            headerString += "21q_cf\t"
            headerString += "22p_cf\t"
            headerString += "22q_cf\t"

            headerString += "1p_cfArmCoverage\t"
            headerString += "1q_cfArmCoverage\t"
            headerString += "2p_cfArmCoverage\t"
            headerString += "2q_cfArmCoverage\t"
            headerString += "3p_cfArmCoverage\t"
            headerString += "3q_cfArmCoverage\t"
            headerString += "4p_cfArmCoverage\t"
            headerString += "4q_cfArmCoverage\t"
            headerString += "5p_cfArmCoverage\t"
            headerString += "5q_cfArmCoverage\t"
            headerString += "6p_cfArmCoverage\t"
            headerString += "6q_cfArmCoverage\t"
            headerString += "7p_cfArmCoverage\t"
            headerString += "7q_cfArmCoverage\t"
            headerString += "8p_cfArmCoverage\t"
            headerString += "8q_cfArmCoverage\t"
            headerString += "9p_cfArmCoverage\t"
            headerString += "9q_cfArmCoverage\t"
            headerString += "10p_cfArmCoverage\t"
            headerString += "10q_cfArmCoverage\t"
            headerString += "11p_cfArmCoverage\t"
            headerString += "11q_cfArmCoverage\t"
            headerString += "12p_cfArmCoverage\t"
            headerString += "12q_cfArmCoverage\t"
            headerString += "13p_cfArmCoverage\t"
            headerString += "13q_cfArmCoverage\t"
            headerString += "14p_cfArmCoverage\t"
            headerString += "14q_cfArmCoverage\t"
            headerString += "15p_cfArmCoverage\t"
            headerString += "15q_cfArmCoverage\t"
            headerString += "16p_cfArmCoverage\t"
            headerString += "16q_cfArmCoverage\t"
            headerString += "17p_cfArmCoverage\t"
            headerString += "17q_cfArmCoverage\t"
            headerString += "18p_cfArmCoverage\t"
            headerString += "18q_cfArmCoverage\t"
            headerString += "19p_cfArmCoverage\t"
            headerString += "19q_cfArmCoverage\t"
            headerString += "20p_cfArmCoverage\t"
            headerString += "20q_cfArmCoverage\t"
            headerString += "21p_cfArmCoverage\t"
            headerString += "21q_cfArmCoverage\t"
            headerString += "22p_cfArmCoverage\t"
            headerString += "22q_cfArmCoverage\t"

            headerString += "1p_cfLevel\t"
            headerString += "1q_cfLevel\t"
            headerString += "2p_cfLevel\t"
            headerString += "2q_cfLevel\t"
            headerString += "3p_cfLevel\t"
            headerString += "3q_cfLevel\t"
            headerString += "4p_cfLevel\t"
            headerString += "4q_cfLevel\t"
            headerString += "5p_cfLevel\t"
            headerString += "5q_cfLevel\t"
            headerString += "6p_cfLevel\t"
            headerString += "6q_cfLevel\t"
            headerString += "7p_cfLevel\t"
            headerString += "7q_cfLevel\t"
            headerString += "8p_cfLevel\t"
            headerString += "8q_cfLevel\t"
            headerString += "9p_cfLevel\t"
            headerString += "9q_cfLevel\t"
            headerString += "10p_cfLevel\t"
            headerString += "10q_cfLevel\t"
            headerString += "11p_cfLevel\t"
            headerString += "11q_cfLevel\t"
            headerString += "12p_cfLevel\t"
            headerString += "12q_cfLevel\t"
            headerString += "13p_cfLevel\t"
            headerString += "13q_cfLevel\t"
            headerString += "14p_cfLevel\t"
            headerString += "14q_cfLevel\t"
            headerString += "15p_cfLevel\t"
            headerString += "15q_cfLevel\t"
            headerString += "16p_cfLevel\t"
            headerString += "16q_cfLevel\t"
            headerString += "17p_cfLevel\t"
            headerString += "17q_cfLevel\t"
            headerString += "18p_cfLevel\t"
            headerString += "18q_cfLevel\t"
            headerString += "19p_cfLevel\t"
            headerString += "19q_cfLevel\t"
            headerString += "20p_cfLevel\t"
            headerString += "20q_cfLevel\t"
            headerString += "21p_cfLevel\t"
            headerString += "21q_cfLevel\t"
            headerString += "22p_cfLevel\t"
            headerString += "22q_cfLevel\t"

            headerString += "largest_cf_level"

            outfile.write(headerString + "\n")

            for curSample in self.runList:
                curLineString = ""
                curLineString += str(curSample.id) + "\t"
                curLineString += str(curSample.cancerType) + "\t"
                curLineString += str(curSample.cancerTypeDetail) + "\t"
                curLineString += str(curSample.purity) + "\t"
                curLineString += str(curSample.onkoCode) + "\t"
                curLineString += str(curSample.tmb) + "\t"
                curLineString += str(curSample.msi) + "\t"
                curLineString += str(curSample.wgd) + "\t"
                curLineString += str(curSample.fga) + "\t"
                curLineString += str(curSample.frac_loh) + "\t"
                curLineString += str(curSample.isHypoploidy) + "\t"
                curLineString += str(curSample.isGainHeavy) + "\t"

                #Add alteration calls.
                for curAlteration in curSample.alterations:
                    curLineString += str(curSample.alterations.get(curAlteration)[0]) + "\t"

                #Add alteration percentages.
                for curAlteration in curSample.alterations:
                    curLineString += str(curSample.alterations.get(curAlteration)[1]) + "\t"

                #Add CF values.
                for i in range(1,23):
                    checkP = str(i) + "p"
                    checkQ = str(i) + "q"
                    pData = curSample.getCFbyArm(checkP)
                    qData = curSample.getCFbyArm(checkQ)
                    curLineString += str(pData[0]) + "\t"
                    curLineString += str(qData[0]) + "\t"

                #Add CF coverage percentages.
                for i in range(1,23):
                    checkP = str(i) + "p"
                    checkQ = str(i) + "q"
                    pData = curSample.getCFbyArm(checkP)
                    qData = curSample.getCFbyArm(checkQ)
                    curLineString += str(pData[1]) + "\t"
                    curLineString += str(qData[1]) + "\t"

                #Get sorted lists of CF data.  
                #Note that this gives us two lists, indexed by chromosome (1 indexed, not 0).
                #One list for p, one list for q.
                sorted_level_data = curSample.getSortedCFs()
                sorted_p = sorted_level_data[0]
                sorted_q = sorted_level_data[1]

                #Add CF level data.
                for i in range(len(sorted_p)):
                    #No such thing as chromosome 0.
                    if i == 0:
                        continue
                    curLineString += str(sorted_p[i]) + "\t"
                    curLineString += str(sorted_q[i]) + "\t"

                curLineString += str(curSample.max_arm_level) + "\t"

                outfile.write(curLineString + "\n")
        
    #This function will print a specific sample in this facets dataset.
    def printFacetsSampleById(self, id):
        try:
            target_sample = self.sampleList[id]
        except KeyError as e:
            print(bcolors.WARNING)
            print(f"\tWarning in FacetsDataset.printFacetsSampleById(): {e}")
            print(bcolors.ENDC)
        else:
            target_sample.printSample()


    #This function will build a facets data set using whatever filters are set.
    def buildFacetsDataset(self):
        if not isinstance(self.ref_facetsMeta, FacetsMeta):
            raise FacetsConfigError("FacetsDataset.buildFacetsDataset(): ref_facetsMeta must be an initialized FacetsMeta.")
        meta = self.ref_facetsMeta

        total_samples_prepped = 0
        num_build_failed      = 0
        num_file_failed       = 0
        num_filtered          = 0

        print("\tApplying filtering and building FacetsDataset beginning with " + str(len(meta.master_file_dict)) + " samples.")

        #Go through our master file dictionary and build facets sample objects.
        for key in list(meta.master_file_dict.keys()):
            if total_samples_prepped % 1000 == 0 and total_samples_prepped != 0:
                print("\t\tSamples Processed: " + str(total_samples_prepped))

            cur_dir_map = meta.master_file_dict.get(key)
            #Clinical maps are keyed by the 17-character tumor sample id, whatever shape the run key has.
            clin_id = short_sample_id(key)

            #CancerType Filter.
            if self.filterCancerType:
                if meta.cancer_type_map.get(clin_id) not in self.selectedCancerTypes:
                    num_filtered += 1
                    continue

            #CancerTypeDetail Filter.
            if self.filterCancerTypeDetail:
                if meta.cancer_type_detail_map.get(clin_id) not in self.selectedCancerDetailTypes:
                    num_filtered += 1
                    continue

            #Clinical Purity Filter.
            if self.filterClinicalPurity:
                curClinPurity = _to_float(meta.clinical_purity_map.get(clin_id))
                if curClinPurity is None or curClinPurity < self.minClinPurity or curClinPurity > self.maxClinPurity:
                    num_filtered += 1
                    continue

            #OnkoCode Filter.
            if self.filterOnkoCode:
                if meta.onkotree_code_map.get(clin_id) not in self.selectedOnkoCodes:
                    num_filtered += 1
                    continue

            #TMB Filter.
            if self.filterTMB:
                cur_tmb = _to_float(meta.cvr_tmb_score_map.get(clin_id))
                if cur_tmb is None or cur_tmb < self.minTMB or cur_tmb > self.maxTMB:
                    num_filtered += 1
                    continue

            #MSI Filter.
            if self.filterMSI:
                cur_msi = _to_float(meta.msi_score_map.get(clin_id))
                if cur_msi is None or cur_msi < self.minMSI or cur_msi > self.maxMSI:
                    num_filtered += 1
                    continue

            #Only parse the .out file once if necessary, not for each individual filter.
            out_data = None
            if self.filterPurity or self.filterPloidy or self.filterDipLogR or self.filterCval or self.filterWGD or self.filterFGA or self.filterFracLoH:
                out_data = FacetsRun.parseOut(cur_dir_map[MetaDictMap.OUT_FILE], meta.hisens_vs_purity)
                if out_data is False:
                    num_file_failed += 1
                    continue

            #Only parse the sampleQC file once if necessary, not for each individual filter.
            sample_qc_data = None
            if self.filterWGD or self.filterFGA or self.filterFracLoH:
                sample_qc_data = FacetsRun.parseSampleQC(cur_dir_map[MetaDictMap.QC_FILE], out_data[1])
                if sample_qc_data is False:
                    num_file_failed += 1
                    continue

            #Purity Filter.
            if self.filterPurity:
                cur_purity = _to_float(out_data[0])
                if cur_purity is None:
                    num_filtered += 1
                    continue
                #FACETS reports purity 0.3 when it could not estimate purity, so that value cannot be trusted.
                if cur_purity == 0.3:
                    num_filtered += 1
                    continue
                if cur_purity < self.minPurity or cur_purity > self.maxPurity:
                    num_filtered += 1
                    continue

            #Ploidy Filter.
            if self.filterPloidy:
                cur_ploidy = _to_float(out_data[2])
                if cur_ploidy is None or cur_ploidy < self.minPloidy or cur_ploidy > self.maxPloidy:
                    num_filtered += 1
                    continue

            #DipLogR Filter.
            if self.filterDipLogR:
                cur_dipLogR = _to_float(out_data[3])
                if cur_dipLogR is None or cur_dipLogR < self.minDipLogR or cur_dipLogR > self.maxDipLogR:
                    num_filtered += 1
                    continue

            #cVal Filter.
            if self.filterCval:
                cur_cVal = _to_float(out_data[1])
                if cur_cVal is None or cur_cVal < self.minC or cur_cVal > self.maxC:
                    num_filtered += 1
                    continue

            #WGD Filter: keepWGD=True keeps only WGD samples, False keeps only non-WGD samples.
            if self.filterWGD:
                if bool(sample_qc_data[0]) != bool(self.keepWGD):
                    num_filtered += 1
                    continue

            #FGA Filter.
            if self.filterFGA:
                if sample_qc_data[1] < self.minFGA or sample_qc_data[1] > self.maxFGA:
                    num_filtered += 1
                    continue

            #Frac LoH Filter.
            if self.filterFracLoH:
                if sample_qc_data[2] < self.minLoh or sample_qc_data[2] > self.maxLoh:
                    num_filtered += 1
                    continue

            #Facets QC Filter: keepQCPass=True keeps only passing fits, False keeps only failing fits.
            if self.filterFacetsQC:
                facets_qc_out = FacetsRun.parseFacetsQC(cur_dir_map)
                if facets_qc_out is None or facets_qc_out == -1:
                    num_file_failed += 1
                    continue
                if bool(facets_qc_out) != bool(self.keepQCPass):
                    num_filtered += 1
                    continue

            #We've passed all filters, we can build the FacetsRun and FacetsSample for this sample.
            try:
                cur_run = FacetsRun.buildFacetsRun(key, meta)
            except ParseError as e:
                meta.failed_samples.append(key)
                if meta.run_verbose:
                    print(bcolors.WARNING + "\t\tWarning: " + str(e) + bcolors.ENDC)
                num_build_failed += 1
                continue
            if cur_run is False:
                num_build_failed += 1
                continue

            #Group runs into FacetsSamples by pair id (<tumor>_<normal>), in both single-run and all-fits mode.
            cur_sample_id = meta.pairIdForKey(key)
            cur_sample = self.sampleList.get(cur_sample_id)
            if cur_sample is None:
                cur_sample = FacetsSample(cur_sample_id)
                self.sampleList[cur_sample_id] = cur_sample
            cur_sample.addRun(cur_run)
            self.runList.append(cur_run)
            total_samples_prepped = total_samples_prepped + 1

        self.numSamples = len(self.runList)
        self.num_file_failed  = num_file_failed
        self.num_build_failed = num_build_failed
        self.num_filtered     = num_filtered
        print("\tFacetsDataset built with " + str(self.numSamples) + " runs in " + str(len(self.sampleList)) + " samples.")
        print("\t\t" + str(num_filtered) + " runs excluded by filters.")
        print("\t\t" + str(num_file_failed) + " runs failed due to file formatting errors.")
        print("\t\t" + str(num_build_failed) + " runs failed during FacetsRun object construction.")


    #This function sets a cancer type filter.  The selectedCancerTypes parameter should be a list of strings.
    def setCancerTypeFilter(self,selectedCancerTypes):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.selectedCancerTypes(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not isinstance(selectedCancerTypes, list):
            raise FacetsConfigError("Error in FacetsDataset.selectedCancerTypes(): selectedCancerTypes must be a list of string values.")
        else:
            self.filterCancerType    = True
            self.selectedCancerTypes = selectedCancerTypes
    #This function sets a cancer type detail filter.  The selectedCancerDetailTypes parameter should be a list of strings.
    def setCancerTypeDetailFilter(self,selectedCancerDetailTypes):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.setCancerTypeDetailFilter(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not isinstance(selectedCancerDetailTypes, list):
            raise FacetsConfigError("Error in FacetsDataset.setCancerTypeDetailFilter(): selectedCancerDetailTypes must be a list of string values.")
        else:
            self.filterCancerTypeDetail    = True
            self.selectedCancerDetailTypes = selectedCancerDetailTypes
    #This function sets a facets purity filter.  
    def setPurityFilter(self, minPurity, maxPurity=1.0):
        if not (isinstance(minPurity, float) or isinstance(minPurity, int)):
            raise FacetsConfigError("Error in FacetsDataset.setCancerTypeDetailFilter(): minPurity must be a numeric value.")
        elif not (isinstance(maxPurity, float) or isinstance(maxPurity, int)):
            raise FacetsConfigError("Error in FacetsDataset.setCancerTypeDetailFilter(): maxPurity must be a numeric value.")
        else:
            self.filterPurity = True
            self.minPurity    = minPurity
            self.maxPurity    = maxPurity
    #This function sets a clinical purity filter.  The selectedCancerDetailTypes parameter should be a list of strings.
    def setClinicalPurityFilter(self, minClinPurity, maxClinPurity=1.0):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.setClinicalPurityFilter(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not (isinstance(minClinPurity, float) or isinstance(minClinPurity, int)):
            raise FacetsConfigError("Error in FacetsDataset.setClinicalPurityFilter(): minClinPurity must be a numeric value.")
        elif not (isinstance(maxClinPurity, float) or isinstance(maxClinPurity, int)):
            raise FacetsConfigError("Error in FacetsDataset.setClinicalPurityFilter(): maxPurity must be a numeric value.")
        else:
            self.filterClinicalPurity = True
            self.minClinPurity        = minClinPurity
            self.maxClinPurity        = maxClinPurity
    #This function sets a onkoCode filter.  The selectedOnkoCodes parameter should be a list of strings.
    def setOnkoCodeFilter(self, selectedOnkoCodes):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.setOnkoCodeFilter(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not isinstance(selectedOnkoCodes, list):
            raise FacetsConfigError("Error in FacetsDataset.setOnkoCodeFilter(): selectedOnkoCodes must be a list of string values.")
        else:
            self.filterOnkoCode    = True
            self.selectedOnkoCodes = selectedOnkoCodes
    #This function sets a ploidy filter.  The min and max ploidy parameters should be a numeric.
    def setPloidyFilter(self, minPloidy, maxPloidy):
        if not (isinstance(minPloidy, float) or isinstance(minPloidy, int)):
            raise FacetsConfigError("Error in FacetsDataset.setPloidyFilter(): minPloidy must be a numeric value.")
        elif not (isinstance(maxPloidy, float) or isinstance(maxPloidy, int)):
            raise FacetsConfigError("Error in FacetsDataset.setPloidyFilter(): maxPloidy must be a numeric value.")
        else:
            self.filterPloidy = True
            self.minPloidy    = minPloidy
            self.maxPloidy    = maxPloidy
    #This function sets a diplogR filter.  The min and max parameters should be a numeric.
    def setDipLogRFilter(self, minDipLogR, maxDipLogR):
        if not (isinstance(minDipLogR, float) or isinstance(minDipLogR, int)):
            raise FacetsConfigError("Error in FacetsDataset.setDipLogRFilter(): minDipLogR must be a numeric value.")
        elif not (isinstance(maxDipLogR, float) or isinstance(maxDipLogR, int)):
            raise FacetsConfigError("Error in FacetsDataset.setDipLogRFilter(): maxDipLogR must be a numeric value.")
        else:
            self.filterDipLogR = True
            self.minDipLogR    = minDipLogR
            self.maxDipLogR    = maxDipLogR
    #This function sets a cVal filter.  The min and max parameters should be a numeric.
    def setCvalFilter(self, minC, maxC):
        if not (isinstance(minC, float) or isinstance(minC, int)):
            raise FacetsConfigError("Error in FacetsDataset.setCvalFilter(): minC must be a numeric value.")
        elif not (isinstance(maxC, float) or isinstance(maxC, int)):
            raise FacetsConfigError("Error in FacetsDataset.setCvalFilter(): maxC must be a numeric value.")
        else:
            self.filterCval = True
            self.minC       = minC
            self.maxC       = maxC
    #This function sets a onkoCode filter.  The keepWGD parameter should be a boolean
    def setWGDFilter(self, keepWGD):
        if not isinstance(keepWGD, bool):
            raise FacetsConfigError("Error in FacetsDataset.setWGDFilter(): keepWGD must be True or False.")
        else:
            self.filterWGD = True
            self.keepWGD   = keepWGD
    #This function sets a FGA filter.  The min and max parameters should be a numeric.
    def setFGAFilter(self, minFGA, maxFGA=1.0):
        if not (isinstance(minFGA, float) or isinstance(minFGA, int)):
            raise FacetsConfigError("Error in FacetsDataset.setFGAFilter(): minFGA must be a numeric value.")
        elif not (isinstance(maxFGA, float) or isinstance(maxFGA, int)):
            raise FacetsConfigError("Error in FacetsDataset.setFGAFilter(): maxFGA must be a numeric value.")
        else:
            self.filterFGA = True
            self.minFGA    = minFGA
            self.maxFGA    = maxFGA
    #This function sets a FGA filter.  The min and max parameters should be a numeric.
    def setFracLohFilter(self, minLoh, maxLoh=1.0):
        if not (isinstance(minLoh, float) or isinstance(minLoh, int)):
            raise FacetsConfigError("Error in FacetsDataset.setFracLohFilter(): minLoh must be a numeric value.")
        elif not (isinstance(maxLoh, float) or isinstance(maxLoh, int)):
            raise FacetsConfigError("Error in FacetsDataset.setFracLohFilter(): maxLoh must be a numeric value.")
        else:
            self.filterFracLoH = True
            self.minLoh        = minLoh
            self.maxLoh        = maxLoh
    #This function sets a facetsQC filter.  The keepQCPass parameter should be a boolean
    def setFacetsQCFilter(self, keepQCPass):
        if not isinstance(keepQCPass, bool):
            raise FacetsConfigError("Error in FacetsDataset.setFacetsQCFilter(): keepQCPass must be True or False.")
        else:
            self.filterFacetsQC = True
            self.keepQCPass     = keepQCPass
    #This function sets a TMB filter.  The min and max parameters should be numeric.
    def setTMBFilter(self, minTMB, maxTMB):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.setTMBFilter(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not (isinstance(minTMB, float) or isinstance(minTMB, int)):
            raise FacetsConfigError("Error in FacetsDataset.setTMBFilter(): minTMB must be a numeric value.")
        elif not (isinstance(maxTMB, float) or isinstance(maxTMB, int)):
            raise FacetsConfigError("Error in FacetsDataset.setTMBFilter(): maxTMB must be a numeric value.")
        else:
            self.filterTMB = True
            self.minTMB    = minTMB
            self.maxTMB    = maxTMB
    #This function sets a MSI filter.  The min and max parameters should be numeric.
    def setMSIFilter(self, minMSI, maxMSI):
        if not self.ref_facetsMeta.hasClinicalMeta:
            raise FacetsConfigError("Error in FacetsDataset.setMSIFilter(): Clinical data was not provided for this run. Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
        if not (isinstance(minMSI, float) or isinstance(minMSI, int)):
            raise FacetsConfigError("Error in FacetsDataset.setMSIFilter(): minMSI must be a numeric value.")
        elif not (isinstance(maxMSI, float) or isinstance(maxMSI, int)):
            raise FacetsConfigError("Error in FacetsDataset.setMSIFilter(): maxMSI must be a numeric value.")
        else:
            self.filterMSI = True
            self.minMSI    = minMSI
            self.maxMSI    = maxMSI
######################
# FacetsGene:    This class represents a single facets gene as represented in the gene_level.txt file.
######################
class FacetsGene:
    def __init__(self, gene, chrom, gene_start, gene_end, seg_start, seg_end, seg_length, cf, tcn, lcn, cn_state, filter, tsg, seg, median_cnlr_seg, segclust, mcn, genes_on_seg, gene_snps, gene_het_snps, spans_segs):
        self.gene            = gene
        self.chrom           = chrom
        self.gene_start      = gene_start
        self.gene_end        = gene_end
        self.seg_start       = seg_start
        self.seg_end         = seg_end
        self.seg_length      = seg_length
        self.cf              = cf
        self.tcn             = tcn
        self.lcn             = lcn
        self.cn_state        = cn_state
        self.filter          = filter
        self.tsg             = tsg
        self.seg             = seg
        self.median_cnlr_seg = median_cnlr_seg
        self.segclust        = segclust
        self.mcn             = mcn
        self.genes_on_seg    = genes_on_seg
        self.gene_snps       = gene_snps
        self.gene_het_snps   = gene_het_snps
        self.spans_segs      = spans_segs



    #Print data about this gene.
    def printGene(self):
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Gene: " + str(self.gene))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Gene Start: " + str(self.gene_start))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Gene End: " + str(self.gene_end))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Seg Start: " + str(self.seg_start))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Seg End: " + str(self.seg_end))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Seg Length: " + str(self.seg_length))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " CF: " + str(self.cf))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " TCN: " + str(self.tcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " LCN: " + str(self.lcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Filter: " + str(self.filter))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " TSG: " + str(self.cf))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Segment: " + str(self.seg))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Median cnlr Segment: " + str(self.median_cnlr_seg))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Segment Cluster: " + str(self.segclust))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " MCN: " + str(self.mcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Genes on Segement: " + str(self.mcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Gene SNPs " + str(self.gene_snps))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Gene het SNPs: " + str(self.gene_het_snps))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Span Segments: " + str(self.spans_segs))
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)

    #Compare two genes.  If they are the same return true.
    def compareGenes(self, geneToCompare):
        if not isinstance(geneToCompare, FacetsGene):
            print (bcolors.FAIL)
            print("Error: FacetsGene.compareSegments(): Can only compare objects of type FacetsGene.")
            print (bcolors.ENDC)
        else:
            if self.gene == geneToCompare.gene and self.gene_start == geneToCompare.gene_start and self.gene_end == geneToCompare.gene_end:
                return True
            else:
                return False
######################
# FacetsSegment:    This class represents a partial segment of facets data as represented in
#                   a facets *.cncf file. cf, tcn, and lcn are all .em values.  
######################
class FacetsSegment:
    percent_arm_loh      = 50.0 # This is the percentage of the arm that needs to be LCN=0 to make an arm level LOH call.
    min_gain_tcn         = 4    # This is the minimum value of TCN required to consider a segment to be a Gain call.
    percent_arm_gain     = 50.0 # This is the percentage of the arm that needs to be TCN=min_gain_tcn to make an arm level Gain call.
    
    def __init__(self, chrom, start, end, cnlr_median, cf, cf_base, tcn, lcn, num_mark, nhet, mafR, segclust, cnlr_median_clust, mafR_clust, adj_mean):
        # These values are calculated in FacetsRun.defineArms().
        self.arm         = "undefined" # This is the chr/arm.  I.E. 1p, 5q, etc.
        self.percentArm  = -1          # This is the percentage of the arm that this segment covers. 

        self.chrom       = int(chrom)
        self.start       = int(start)
        self.end         = int(end)
        self.cnlr_median = float(cnlr_median) #These values are the same as the unadjusted.seg file so we did not use the unadjusted.seg file to open less files
        self.length      = int(max(end - start, start - end))
        self.cf          = float(cf)          #This API uses cf.em values for everything, so its just called cf here.
        self.num_mark    = int(num_mark)
        self.nhet        = int(nhet)
        self.segclust    = int(segclust)
        self.cnlr_median_clust= float(cnlr_median_clust) 
        self.mafR  = float(mafR)
        self.adj_mean = float(adj_mean) #These values come from the adjusted.seg file

        try:
            self.cf_base = float(cf_base)
        except ValueError:
            try:
                self.cf_base = int(cf_base)
            except ValueError:
                self.cf_base = -1

        if pd.isna(mafR_clust):
            self.mafR_clust     = 0
        else:
            self.mafR_clust     = float(mafR_clust)


        if pd.isna(tcn):
            self.tcn     = -1
        else:
            self.tcn     = int(tcn)
        if pd.isna(lcn):
            self.lcn     = -1
        else:
            self.lcn     = int(lcn)
        self.isLoH       = False
        self.isGain      = False

    #Compare two segments.  If they are the same return true.
    def compareSegments(self, segToCompare):
        if not isinstance(segToCompare, FacetsSegment):
            print (bcolors.FAIL)
            print("Error: Can only compare objects of type FacetsSegment.")
            print (bcolors.ENDC)
        else:
            if self.chrom == segToCompare.chrom and self.start == segToCompare.start and self.end == segToCompare.end:
                return True
            else:
                return False
    #Print a summary of this segment.
    def printSegment(self):
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Chromosome: " + str(self.chrom))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Arm: " + (self.arm))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Percent Arm: " + str(self.percentArm))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Start Position: " + str(self.start))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " End Position: " + str(self.end))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Length: " + str(self.length))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " CF (base): " + str(self.cf_base))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " CF (em): " + str(self.cf))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " TCN (em): " + str(self.tcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " LCN (em): " + str(self.lcn))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Number Mark: " + str(self.num_mark))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Number Het: " + str(self.nhet))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Unadj Median: " + str(self.cnlr_median)) 
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Adj mean: " + str(self.adj_mean))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " mafR: " + str(self.mafR))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " Segment Clust: " + str(self.segclust))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " cnlr_median_clust: " + str(self.cnlr_median_clust))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " isLoH: " + str(self.isLoH))
        print(bcolors.BOLD + "\t|" + bcolors.ENDC + " isGain: " + str(self.isGain))
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)


######################
# FacetsSample:    This class represents a FacetsSample containing one or more runs. 
######################
class FacetsSample:
    def __init__(self, id):
        self.id   = id
        self.runs = []
    
    #Print information about this sample.
    def printSample(self):
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)
        print("FacetsSample: " + str(self.id))
        print("\tThis sample contains " + str(len(self.runs)) + " FacetsRuns.")
        for i in range(len(self.runs)):
            print("ID: " + str(self.runs[i].id))
            print("Fit Dir: " + str(self.runs[i].fitDir))
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)

    #Add a FacetsRun to this FacetsSample.
    def addRun(self, runToAdd):
        if not isinstance(runToAdd, FacetsRun):
            raise FacetsConfigError("Error: FacetsSample.addRun(): Can only add objects of type FacetsRun.")
        else:
            self.runs.append(runToAdd)
            return True

######################
# FacetsRun:     This class represents single run for a facets sample, and contains a variety of 
#                   metadata for the run, as well as a list of FacetsSegments and FacetsGenes.  
######################
class FacetsRun:
    alteration_arm_num = 20    # The number of arms that need to be LOH or Gain to make a sample level Hypoploidy or GainHeavy call.
    alteration_cov_min = 50    # The percentage of the arm that needs to be covered to make LOH/Gain calls.
    run_alterations    = False # Whether or not to run the alteration calculation functions.
    split_arms         = False # Wheter or not to split segments into p and q arms.

    def __init__(self, id, patientId, fitDir, cancerType, cancerTypeDetail, 
                 purity, clinicalPurity, onkoCode, ploidy, dipLogR, cval, 
                 wgd, fga, frac_loh, facets_qc, tmb, msi, purity_hisens,
                 review_status, fit_name, review_author, fit_class="standard"):
        self.segments         = []
        self.genes            = []
        self.id               = id
        self.fit_class        = fit_class   # "standard", or "clinical"/"research" for a 2n pair
        self.patientId        = patientId
        self.fitDir           = fitDir
        self.cancerType       = cancerType
        self.cancerTypeDetail = cancerTypeDetail
        self.purity           = purity
        self.clinicalPurity   = clinicalPurity
        self.onkoCode         = onkoCode
        self.ploidy           = ploidy
        self.dipLogR          = dipLogR
        self.cval             = cval
        self.wgd              = bool(wgd)
        self.fga              = fga
        self.frac_loh         = frac_loh
        self.facets_qc        = facets_qc
        self.tmb              = tmb
        self.msi              = msi
        self.purity_hisens    = purity_hisens
        self.review_status    = review_status
        self.fit_name         = fit_name
        self.review_author    = review_author
        self.alterations      = {} #A map of {chr_arm -> [alteration, percentAltered]}, where alteration is (LoH or Gain or Neutral) based on the segments for this sample.
        self.isHypoploidy     = False
        self.isGainHeavy      = False
        self.cfLevel_arms     = {} #A map of {cf_level -> [chr_arms]}
        self.cfLevel_values   = {} #A map of {cf_level -> cf_value}
        self.arm_cf_data      = {} #A map of {arm -> [cf, percent_arm_coverage]}
        self.max_arm_level    = -1 #The arm level that has the most members.
        

    #Get CN_call value for any given tcn/lcn.
    def get_cn_call(self, tcn, lcn):
            """
            Return the copy-number call string based on this run's WGD status,
            total copy number (tcn) and lesser copy number (lcn).
            """
            # NA handling (None / NaN / the -1 sentinel) and the lookup live in cn_calls.get_cn_call.
            return _get_cn_call(self.wgd, tcn, lcn)

    #Add a segment to the segment array of type FacetsSegment.
    def addSegment(self, seg):
        if not isinstance(seg, FacetsSegment):
            raise FacetsConfigError("Error: FacetsRun.addSegment(): Can only add objects of type FacetsSegment.")
        else:
            self.segments.append(seg)
            return True

    #Remove a FacetsSegment from the segment array. Returns True if successful.
    def removeSegment(self, segToRemove):
        startingNumSegs = len(self.segments)
        for i in range(len(self.segments)):
            if self.segments[i].compareSegments(segToRemove):
                self.segments.pop(i)
                if len(self.segments) == startingNumSegs - 1:
                    return True
                else:
                    raise FacetsAPIError("Error in FacetsRun.removeSegment(): Remove failed to shrink segment array.")
        return False

    #Print all genes in the gene array.
    def printAllGenes(self):
        if len(self.genes) == 0:
            print (bcolors.WARNING)
            print("Warning in FacetsRun.printSegments(): No segments available for sample " + self.id)
            print (bcolors.ENDC)
        else:
            for i in range(len(self.segments)):
                print(bcolors.BOLD + "\t|-----------| Gene " + str(i) + " |-----------|" + bcolors.ENDC)
                self.genes[i].printGene()

    #Add a FacetsGene to the genes array.  Returns True if successful.
    def addGene(self,gene):
        if not isinstance(gene, FacetsGene):
            raise FacetsConfigError("Error: FacetsRun.addGene(): Can only add objects of type FacetsGene.")
        else:
            self.genes.append(gene)
            return True

    #Remove a FacetsGene from the genes array. Returns True if successful.
    def removeGene(self, geneToRemove):
        startingNumGenes = len(self.genes)
        for i in range(len(self.genes)):
            if self.genes[i].compareGenes(geneToRemove):
                self.genes.pop(i)
                if len(self.genes) == startingNumGenes - 1:
                    return True
                else:
                    raise FacetsAPIError("Error in FacetsRun.removeGene(): Remove failed to shrink gene array.")
        return False

    #This function, given a target chromosome and target arm, will return a list of all segments for this sample that fall in the requested regions.
    def getSegmentsByChromAndArm(self, target_chrom, target_arm):
        if(int(target_chrom) <= 0 or int(target_chrom) >= 23):
            raise FacetsConfigError("Error in FacetsRun.getSegmentsByChromAndArm(): target_chrom must be between 1 and 23.")
        if(target_arm != "p" and target_arm != "q"):
            raise FacetsConfigError("Error in FacetsRun.getSegmentsByChromAndArm(): target_arm must be 'p' or 'q'.")

        arm_to_fetch = str(target_chrom) + str(target_arm)
        
        target_seg_list = []
        for i in range(len(self.segments)):
            #self.segments[i].printSegment()
            if self.segments[i].arm == arm_to_fetch:
                target_seg_list.append(self.segments[i])

        return target_seg_list 
        
    #This function, given a target chromosome, will return a list of all segments for this sample that fall in the requested regions.
    def getSegmentsByChrom(self, target_chrom):
        if(target_chrom <= 0 or target_chrom >= 23):
            raise FacetsConfigError("Error in FacetsRun.getSegmentsByChrom(): target_chrom must be between 1 and 23.")

        target_seg_list = []
        for i in range(len(self.segments)):
            if self.segments[i].chrom == target_chrom:
                target_seg_list.append(self.segments[i])

        return target_seg_list 
        
    #This function will calculate purity based on all of the segments associated with this sample.
    #Purity is the maximum cf value that is not 1 or NA. if use_base is set to true, this function
    #will calculate the purity based on the non .em cf values.  Otherwise cf.em will be used.
    def calculatePurityByCF(self, use_base=False):
        all_cfs      = []
        purity_value = -1
        for cur_seg in self.segments:
            if use_base:
                if cur_seg.cf_base == -1 or cur_seg.cf_base == 1:
                    continue
                all_cfs.append(cur_seg.cf_base)
            else:
                if cur_seg.cf == -1 or cur_seg.cf == 1:
                    continue
                all_cfs.append(cur_seg.cf)
        if(len(all_cfs) == 0):
            return 1
        else:
            purity_value = max(all_cfs)
        
        return purity_value

    #This function will return the number of arms with LOH over 'percent_lost' in this sample.
    def getNumberLossArms(self, percent_lost):
        total_lost = 0
        for arm in self.alterations:
            #Alterations: {chr_arm -> [alteration, percentAltered]}
            cur_arm = self.alterations.get(arm)
            if (cur_arm[0] == "LOH" or cur_arm[0] == "LOH/GAIN") and cur_arm[1] >= percent_lost:
                total_lost = total_lost + 1

        return total_lost

    #This function will return the number of arms with GAIN over 'percent_gain' in this sample.
    def getNumberGainArms(self, percent_gain):
        total_gain = 0
        for arm in self.alterations:
            #Alterations: {chr_arm -> [alteration, percentAltered]}
            cur_arm = self.alterations.get(arm)
            if (cur_arm[0] == "GAIN" or cur_arm[0] == "LOH/GAIN") and cur_arm[1] >= percent_gain:
                total_gain = total_gain + 1

        return total_gain

    #This function will sort the self.cfLevel_arms into and return an indexed list of arm levels.
    #The function returns two lists, 1 indexed (ignore 0 position), for [[p_arms],[q_arms]].
    def getSortedCFs(self):
        if len(self.cfLevel_arms) == 0:
            raise FacetsAPIError("FacetsRun.getSortedCFs(): cfLevels is empty. Run defineCFLevels() first.")

        #Keep two lists: one for p, one for q. 
        #Populate them initially with -1s so we keep track of anything missing.
        #Make it one longer so we can 1 index this instead of 0 index.
        p_list = []
        q_list = []
        for i in range(0,23):
            p_list.append(-1)
            q_list.append(-1)

        #Assign arms to their correct location in our lists.
        for curLevel in self.cfLevel_arms:
            armsInLevel = self.cfLevel_arms.get(curLevel)
            for curChrArm in armsInLevel:
                curChr = int(curChrArm[:-1])
                curArm = curChrArm[-1]
                if curArm == "p":
                    p_list[curChr] = curLevel
                if curArm == "q":
                    q_list[curChr] = curLevel

        return [p_list, q_list]

    #This function returns cf data for a specified arm, i.e. 1p.
    #It is used because arm_cf_data isn't sorted since its a dictionary and sometimes we want to work
    #with data in a specific order, i.e. when writing files.
    def getCFbyArm(self, target_arm):
        return self.arm_cf_data.get(target_arm)
    # This function will calculate a CF value covering the largest segment of a target_arm.
    # I.E. if an arm is 70% cf=.52 and 30% cf=.23, this function will return .52.  
    # match_cf_range is the maximum distance that two segments can be from one another in order to be considered the same CF.
    # This function returns [longest_cf_value, percent_of_arm_covered].
    def getArmLevelCF(self, target_arm, match_cf_range=0.01):
        if target_arm == "undefined":
            return False

        curArmSegs = self.getSegmentsByChromAndArm(target_arm[:-1], target_arm[-1])
        #If there are no segments for this arm, we want to return some negative value instead of just crashing.
        if len(curArmSegs) == 0:
            return False

        cf_percents = {} # A map of cf values and how much they cover. i.e: .37 -> .25, .66 -> .65
        for curSeg in curArmSegs:
            if curSeg.arm == "undefined" or curSeg.percentArm == -1:
                raise FacetsAPIError(str("Error in FacetsRun.getArmLevelCF(): Undefined arms in " + self.id + ". Run defineArms()."))

            matched_seg = False
            for observed_cf in cf_percents:
                cf_distance = abs(observed_cf - curSeg.cf)
                #If we find a segment that has a matching CF with another segment, add to that CF's total coverage percent.
                if cf_distance <= match_cf_range:
                    cur_cf_percent = cf_percents.get(observed_cf)
                    cf_percents[observed_cf] = cur_cf_percent + curSeg.percentArm
                    matched_seg = True
            #If there are no other CFs matching our current segment, add a new CF to the cf_percents.
            if not matched_seg:
                cf_percents[curSeg.cf] = curSeg.percentArm

        #Now find the largest percent coverage
        max_cf = -1.0
        arm_coverage_percent = -1.0
        for observed_cf in cf_percents:
            cur_cf_percent = cf_percents.get(observed_cf)
            if arm_coverage_percent < cur_cf_percent:
                max_cf = observed_cf
                arm_coverage_percent = cur_cf_percent
            
        if max_cf == -1:
            return False  # We didn't find any arms with a large enough CF to make an arm call.
        else:
            return [max_cf, arm_coverage_percent] # We return the cf that has the most coverage for this arm and what percent coverage it had.

    #This function is for merging CF levels.  Helper for defineCFLevels(). 
    #Given a list of levels to merge [[1,3], [2,4], [2,3]] this function should merge
    #the sub-lists into sets where any overlap between sets is given.  In the above example this should return [[1,2,3,4]].
    def mergeCFLevels(self, lsts):
        sets = [set(lst) for lst in lsts if lst]
        merged = True
        while merged:
            merged = False
            results = []
            while sets:
                common, rest = sets[0], sets[1:]
                sets = []
                for x in rest:
                    if x.isdisjoint(common):
                        sets.append(x)
                    else:
                        merged = True
                        common |= x
                results.append(common)
            sets = results
        return sets

    # Use defined arms and determine the cf levels across all arms.  
    # A CF level is a CF value that matches other CF levels in a sample. All samples with the same CFs will be in the same level.
    # match_cf_range is the maximum distance that two segments can be from one another in order to be called at the same level.
    # This function will identify arms that belong on the same level together and assign data for CF values and maps to this
    # FacetsRun object.
    def defineCFLevels(self, match_cf_range=0.01):
        if len(self.segments) == 0:
            raise FacetsAPIError(str("Error in FacetsRun.defineCFLevels(): No segments available for sample " + self.id))
        if len(self.alterations) == 0:
            raise FacetsAPIError(str("Error in FacetsRun.defineCFLevels(): Alterations have not been defined for " + self.id + ". Run defineArms() first."))
        
        existingLevels    = {}  # This is to keep track of how many levels already exist. I.E: {1 -> ["1p","3q"], 2 -> ["5p", "5q", "12p"]}
        levelCFs          = {}  # This is to keep track of the CF value for each given level. {1 -> .35, 2 -> .27}
        curMaxLevel       = 1   # This is the highest existing level at any given time.  Increment as we add more levels.
        
        #Make a list of which arms we need to check.
        arms_to_check = []   # The arms that we have segment data for in this sample.
        for curSeg in self.segments:
            if curSeg.arm not in arms_to_check:
                arms_to_check.append(curSeg.arm)

        #Get CF data for all of our arms.  {1p -> [cur_arm_cf, percent_arm_covered], [...]}
        all_arm_cfs = {} 
        for curArm in arms_to_check:             
            if self.getArmLevelCF(curArm) is False:
                continue
            all_arm_cfs[curArm] = self.getArmLevelCF(curArm)

        #Find all of the unique CF values.
        initial_levels = []
        for base_arm in all_arm_cfs:
            base_arm_cf = all_arm_cfs.get(base_arm)[0]
            if base_arm_cf not in initial_levels:
                initial_levels.append(base_arm_cf)

        #Identify pairwise sets of elements from our initial levels that should be merged based on close/identical cf distance.
        levels_to_merge = []
        for i in range(len(initial_levels)):
            for j in range(len(initial_levels)):
                if i == j:
                    continue
                cfDist = abs(initial_levels[i] - initial_levels[j])
                if cfDist <= match_cf_range:
                    levels_to_merge.append([i,j])

        #Merge the pairwise sets into unique level sets.
        merged_levels = self.mergeCFLevels(levels_to_merge)

        #For our final map, we don't want to have duplicates.  
        #We can just look at the first element in any level and know they share a CF at that point.
        #Build up levels_to_ignore to know which arms we dont need to look at once we have a CF sorted out for that level.
        levels_to_ignore = []
        for i in range(len(merged_levels)):
            for j in range(len(merged_levels[i])):
                if j == 0:
                    continue
                else:
                    levels_to_ignore.append(list(merged_levels[i])[j])

        #Build the levels.  We want to have a set of {level -> cf_value} so we know which level has which CF.
        #All possible CFs for the sample should be assigned to a level within a range specified by match_cf_range.
        for i in range(len(initial_levels)):
            if i in levels_to_ignore:
                continue
            levelCFs[curMaxLevel] = initial_levels[i]
            curMaxLevel = curMaxLevel + 1
        
        #Now we can look over all of our arm level cf data, and compare it with our levelCFs map.
        #If the observed level for a given sample matches a specific level (which it must at this point),
        #we add that arm to the specific level for a map of {level -> [arm_X,arm_Y]}
        for curLevel in levelCFs:
            curLevel_cf = levelCFs.get(curLevel)
            for curArm in all_arm_cfs:
                curArm_cf = all_arm_cfs.get(curArm)[0]
                
                cfDist = abs(curLevel_cf - curArm_cf)
                if cfDist <= match_cf_range:
                    if curLevel not in existingLevels:
                        existingLevels[curLevel] = [curArm]
                    else:
                        curExistingLevel = existingLevels.get(curLevel)
                        curExistingLevel.append(curArm)

        #Clean up the all_arm_cfs so that it has values for every chromosome, 
        #even if there were no segments present for a given arm.
        armsToAdd = []
        for i in range(1,23):
            checkP = str(i) + "p"
            checkQ = str(i) + "q"
            if checkP not in all_arm_cfs:
                armsToAdd.append(checkP)
            if checkQ not in all_arm_cfs:
                armsToAdd.append(checkQ) 

        for item in armsToAdd:
            all_arm_cfs[item] = [-1,-1]
            

        #We need to check now that every chromosome is present in the existing levels.
        #Sometimes it seems like singleton levels get ignored and marked -1, so this section adjusts that.
        allExistingLevels = []
        for curLevel in existingLevels:
            levelVals = existingLevels.get(curLevel)
            for item in levelVals:
                allExistingLevels.append(item)

        missingLevels = []
        for i in range(1,23):
            checkP = str(i) + "p"
            checkQ = str(i) + "q"
            if checkP not in allExistingLevels:
                missingLevels.append(checkP)
            if checkQ not in allExistingLevels:
                missingLevels.append(checkQ) 

        for item in missingLevels:
            if item in all_arm_cfs:
                if all_arm_cfs.get(item)[0] == -1:
                    continue
                existingLevels[curMaxLevel] = [item]
                curMaxLevel = curMaxLevel + 1
        
        #Identify what the largest CF level is.
        largest_level = -1
        largest_level_len = -1
        for curLevel in existingLevels:
            if len(existingLevels.get(curLevel)) > largest_level_len:
                largest_level_len = len(existingLevels.get(curLevel))
                largest_level = curLevel

        #Now we can add all of these data to their corresponding FacetsRun level places for future processing.
        self.cfLevel_arms   = existingLevels
        self.cfLevel_values = levelCFs
        self.arm_cf_data    = all_arm_cfs
        self.max_arm_level  = largest_level

    #Check segments to determine p or q arm. If centromere is crossed, split the segment into two new segments.
    def defineArms(self):
        if len(self.segments) == 0:
            print (bcolors.WARNING)
            print("Warning in FacetsRun.defineArms(): No segments available for sample " + self.id)
            print (bcolors.ENDC)
        else:
            segsToRemove = []
            segsToAdd    = []
            for i in range(len(self.segments)):
                cur_seg = self.segments[i]

                #We don't consider chromosome 23.
                if cur_seg.chrom == 23:
                    continue

                if cur_seg.start > cur_seg.end:
                    raise ParseError("FacetsRun.defineArms(): end position occurs before start position in %s (chrom %s: %s-%s)."
                                     % (self.id, cur_seg.chrom, cur_seg.start, cur_seg.end))

                #Determine if the centromere is intersected.
                cur_chrom_ranges = FacetsMeta.chr_arms.get(cur_seg.chrom)
                lowPos  = min(cur_seg.start, cur_seg.end)
                highPos = max(cur_seg.start, cur_seg.end)
                centromere_pos = cur_chrom_ranges[2]

                #If everything occurs in the P or q arm range only, assign the arm.
                if lowPos <= centromere_pos and highPos <= centromere_pos:
                    cur_seg.arm = str(cur_seg.chrom) + "p"
                    cur_seg.percentArm = float(cur_seg.length) / float(cur_chrom_ranges[1])
                elif lowPos >= centromere_pos and highPos >= centromere_pos:
                    cur_seg.arm = str(cur_seg.chrom) + "q"
                    cur_seg.percentArm = float(cur_seg.length) / (float(cur_chrom_ranges[3]) - float(cur_chrom_ranges[2]))
                #If not, split the segment.
                else:
                    segsToRemove.append(cur_seg)
                    
                    #Split into two segments around the centromere.
                    seg1_newStart = cur_seg.start
                    seg1_newEnd   = cur_chrom_ranges[1] - 1
                    new_Seg1      = FacetsSegment(cur_seg.chrom,seg1_newStart,seg1_newEnd,cur_seg.cnlr_median,cur_seg.cf,cur_seg.cf_base,cur_seg.tcn,cur_seg.lcn,cur_seg.num_mark, cur_seg.nhet, cur_seg.mafR, cur_seg.segclust, cur_seg.cnlr_median_clust, cur_seg.mafR_clust, cur_seg.adj_mean)
                    new_Seg1.arm  = str(cur_seg.chrom) + "p"
                    new_Seg1.percentArm = float(new_Seg1.length) / float(cur_chrom_ranges[1])

                    seg2_newStart = cur_chrom_ranges[2] + 1
                    seg2_newEnd   = cur_seg.end
                    new_Seg2      = FacetsSegment(cur_seg.chrom,seg2_newStart,seg2_newEnd,cur_seg.cnlr_median,cur_seg.cf,cur_seg.cf_base,cur_seg.tcn,cur_seg.lcn,cur_seg.num_mark, cur_seg.nhet, cur_seg.mafR, cur_seg.segclust, cur_seg.cnlr_median_clust, cur_seg.mafR_clust, cur_seg.adj_mean)
                    new_Seg2.arm  = str(cur_seg.chrom) + "q"
                    new_Seg2.percentArm = float(new_Seg2.length) / (float(cur_chrom_ranges[3]) - float(cur_chrom_ranges[2]))

                    segsToAdd.append(new_Seg1)
                    segsToAdd.append(new_Seg2)

            #Remove segments that needed to be split.
            for i in range(len(segsToRemove)):
                removeOk = self.removeSegment(segsToRemove[i])
                if not removeOk:
                    raise FacetsAPIError("Error in FacetsRun.defineArms(). Failed to remove centromere split segment.")

            #Add segments created from splitting old segments.
            for i in range(len(segsToAdd)):
                addOk = self.addSegment(segsToAdd[i])
                if not addOk:
                    raise FacetsAPIError("Error in FacetsRun.defineArms(). Failed to add centromere split segment.")
        
    #This function will check this sample and determine if it is a hypoploidy sample.
    #It will make these determinations based on a minimum percent of arm covered for min_arm_num of samples.
    #This function will set isHypoploidy and isGainHeavy for this sample.
    def defineSampleLevelAlterationCall(self, min_percent_arm_cov, min_arm_num):
        if self.getNumberLossArms(min_percent_arm_cov) >= min_arm_num:
            self.isHypoploidy = True
        if self.getNumberGainArms(min_percent_arm_cov) >= min_arm_num:
            self.isGainHeavy = True
    #This is a helper function for defineAllLohAndGains().  
    #It will perform LoH and Gain calculations on a set of segments for a given chromosome and arm.
    def calculateLohAndGains(self, selected_segs, target_chr, target_arm):
        total_lcn0_len        = 0
        total_tcnGain_len     = 0
        total_lcn0_percent    = 0
        total_tcnGain_percent = 0

        #Adjust Gain calls if this is a WGD sample.
        wgd_adjusted_min_gain_tcn = FacetsSegment.min_gain_tcn
        if self.wgd:
            wgd_adjusted_min_gain_tcn = wgd_adjusted_min_gain_tcn + 2
        
        #Determine length of whatever chromosome arm we are looking at.
        cur_chr_arm           = str(target_chr) + target_arm
        cur_chr_arms          = FacetsMeta.chr_arms.get(target_chr)
        cur_arm_len           = 0
        if target_arm == "p":
            cur_arm_len = cur_chr_arms[1] - cur_chr_arms[0]
        if target_arm == "q":
            cur_arm_len = cur_chr_arms[3] - cur_chr_arms[2]

        for x in range(len(selected_segs)):
            #Don't even look at segments that don't meet our LoH/Gain criteria.
            if selected_segs[x].lcn > 0 and selected_segs[x].tcn < wgd_adjusted_min_gain_tcn:
                continue
            else:
                #For the N/A LCNs, just move on as well.
                if selected_segs[x].lcn == -1:
                    continue

                #Count up the total lengths of segments meeting LoH or Gain criteria.
                if selected_segs[x].lcn <= 0:
                    total_lcn0_len = total_lcn0_len + selected_segs[x].length
                if selected_segs[x].tcn >= wgd_adjusted_min_gain_tcn:
                    total_tcnGain_len = total_tcnGain_len + selected_segs[x].length

        #Calculate LoH and Gain percentages.
        if total_lcn0_len > cur_arm_len or total_tcnGain_len > cur_arm_len:
            raise FacetsAPIError("Error in FacetsRun.calculateLohAndGains(). Total LCN=0 length or TCN>=Gain is greater than arm length.")
        else:
            total_lcn0_percent    = float((total_lcn0_len / cur_arm_len) * 100.0)
            total_tcnGain_percent = float((total_tcnGain_len / cur_arm_len) * 100.0)

            #Each chr/arm combination must be unique in the alterations map.
            if cur_chr_arm in self.alterations:
                raise FacetsAPIError("Error in FacetsRun.calculateLohAndGains(). Duplicate key in alteration map.")

            #Update the segments themselves and the sample's alterations map.
            if total_lcn0_percent >= FacetsSegment.percent_arm_loh:
                selected_segs[x].isLoH = True
                self.alterations[cur_chr_arm] = ["LOH", total_lcn0_percent]
            elif total_tcnGain_percent >= FacetsSegment.percent_arm_gain:
                selected_segs[x].isGain = True
                self.alterations[cur_chr_arm] = ["GAIN", total_tcnGain_percent]
            elif total_lcn0_percent >= FacetsSegment.percent_arm_loh and total_tcnGain_percent >= FacetsSegment.percent_arm_gain:
                selected_segs[x].isGain = True
                selected_segs[x].isLoH = True
                self.alterations[cur_chr_arm] = ["LOH/GAIN", total_lcn0_percent]
            else:
                self.alterations[cur_chr_arm] = ["NEUTRAL", 0]

    ######################
    # defineAllLohAndGains:  This function will process all segments for this FacetsRun object and populate
    #                        the 'alterations' map for the sample.  A sample is determined to be LoH if it is 0 LCN for
    #                        more than the defined FacetsSegment threshold percent of of the arm for all segments on that arm.
    #                        A gain is determined by a correspondingly TCN greater than FacetsSegment.min_gain_tcn for FacetsSegment.percent_arm_gain.
    ######################
    def defineAllLohAndGains(self):
        for i in range(1,23):
            #Handle p and q arms for each chromosome.
            cur_segs_p    = self.getSegmentsByChromAndArm(i, "p")
            cur_segs_q    = self.getSegmentsByChromAndArm(i, "q")
            self.calculateLohAndGains(cur_segs_p, i, "p")
            self.calculateLohAndGains(cur_segs_q, i, "q")

    #Print all segments in the segment array.
    def printAllSegments(self):
        if len(self.segments) == 0:
            print (bcolors.WARNING)
            print("Warning in FacetsRun.printSegments(): No segments available for sample " + self.id)
            print (bcolors.ENDC)
        else:
            for i in range(len(self.segments)):
                print(bcolors.BOLD + "\t|-----------| Segment " + str(i) + " |-----------|" + bcolors.ENDC)
                self.segments[i].printSegment()


    #Print a specific segment by index.
    def printTargetSegment(self, segNumber):
        if len(self.segments) == 0:
            print (bcolors.WARNING)
            print("Warning in FacetsRun.printSegments(): No segments available for sample " + self.id)
            print (bcolors.ENDC)
        else:
            print(bcolors.BOLD + "\t|-----------| Segment " + str(segNumber) + " |-----------|" + bcolors.ENDC)
            self.segments[segNumber].printSegment()


    #Print a summary of this sample.
    def printSample(self):
        print(bcolors.BOLD + "~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Patient: " + self.patientId)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Sample:  " + self.id)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Directory: " + self.fitDir)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Purity/Hisens: " + self.purity_hisens)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Passed Facets QC: " + str(self.facets_qc))
        print(bcolors.BOLD + "| --------------------------------------- ")
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Cancer Type: " + self.cancerType)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Cancer Type Detailed: " + self.cancerTypeDetail)
        print(bcolors.BOLD + "|" + bcolors.ENDC + " OnkoTree Code: " + self.onkoCode)
        print(bcolors.BOLD + "| --------------------------------------- ")
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Clinical Purity: " + str(self.clinicalPurity))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " FACETS Purity: " + str(self.purity))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Ploidy: " + str(self.ploidy))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " cVal: " + str(self.cval))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " DipLogR: " + str(self.dipLogR))
        print(bcolors.BOLD + "| --------------------------------------- ")
        print(bcolors.BOLD + "|" + bcolors.ENDC + " Fraction LoH: " + str(self.frac_loh))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " FGA: " + str(self.fga))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " WGD: " + str(self.wgd))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " TMB: " + str(self.tmb))
        print(bcolors.BOLD + "|" + bcolors.ENDC + " MSI: " + str(self.msi))
        print(bcolors.BOLD + "| --------------------------------------- ")
        print(bcolors.BOLD + "|" + bcolors.ENDC + " This FacetsRun currently has " + str(len(self.segments)) + " associated segments.")
        print(bcolors.BOLD + "~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC )


    ######################
    # parseOut:  This function accepts a path to a .out file from a facets directory.
    #            It will extract and return purity, cval, ploidy, and dipLogR values as a list of strings
    #            (as written in the file; "" when a key is absent). Returns False if the file cannot be read.
    #            mode selects the cval line: "purity" -> '# purity_cval', "hisens" -> '# cval'. When omitted,
    #            the FacetsMeta class default is used (backward compatibility); FacetsMeta always passes its own.
    #            Returns: [purity,cval,ploidy,dipLogR]
    ######################
    @staticmethod
    def parseOut(out_file, mode=None):
        mode = validate_mode(mode if mode is not None else FacetsMeta.hisens_vs_purity)
        cval_key = "# purity_cval" if mode == "purity" else "# cval"
        vals = {"# Purity": "", cval_key: "", "# Ploidy": "", "# dipLogR": ""}
        try:
            with open(out_file, 'r') as fh:
                for line in fh:
                    if "=" not in line:
                        continue
                    key, _sep, value = line.partition("=")
                    key = key.strip()
                    if key in vals and vals[key] == "":
                        vals[key] = value.strip()
        except OSError as e:
            if FacetsMeta.run_verbose:
                print(bcolors.WARNING + "\t\tCannot read .out file %s: %s" % (out_file, e) + bcolors.ENDC)
            return False
        return [vals["# Purity"], vals[cval_key], vals["# Ploidy"], vals["# dipLogR"]]

    ######################
    # getSNPs:  This function accepts a path to a .Rdata file and returns SNP data (out$jointseg) 
    #           as a rpy2 table.
    ######################
    #@staticmethod
    #def getSNPs(rdata_file):
    #    try:
    #        robjects.r['load'](rdata_file)
    #        r_data = robjects.r['out']

            # Convert R object to pandas DataFrame
    #        pandas2ri.activate()
    #        df = pandas2ri.ri2py(r_data.rx2('jointseg'))


    ######################
    # parseGeneLevel:  This function accepts a path to a .gene_level.txt file and returns a list of lists
    #                  containing data representing the gene level file.
    #                  returns [[gene1,start,end...],[gene2,start,end...]]
    ######################
    @staticmethod
    def parseGeneLevel(gene_level_file):
        try:
            #Extract Gene Level data.
            gene_level_df       = pd.read_csv(gene_level_file, sep="\t", low_memory=False)
            cur_gene_level_data = []
            for index in gene_level_df.index:
                cur_gene_level_data.append([gene_level_df['gene'][index], gene_level_df['chrom'][index], gene_level_df['gene_start'][index], gene_level_df
                                            ['gene_end'][index], gene_level_df['seg_start'][index], gene_level_df['seg_end'][index], 
                                            gene_level_df['seg_length'][index], gene_level_df['cf.em'][index], gene_level_df['tcn.em'][index], 
                                            gene_level_df['lcn.em'][index], gene_level_df['cn_state'][index], gene_level_df['filter'][index],
                                            gene_level_df['tsg'][index], gene_level_df['seg'][index], gene_level_df['median_cnlr_seg'][index],
                                            gene_level_df['segclust'][index], gene_level_df['mcn'][index], gene_level_df['genes_on_seg'][index],
                                            gene_level_df['gene_snps'][index], gene_level_df['gene_het_snps'][index], gene_level_df['spans_segs'][index],
                                            ])
            
            return cur_gene_level_data

        except (OSError, KeyError, ValueError, TypeError, IndexError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            raise ParseError("Error parsing sample qc file. Terminating execution." + " " + str("\t" + gene_level_file) + ": " + str(e)) from e
    ######################
    # parseSampleQC:  This function accepts a path to a .qc.txt file and an expected cval
    #                 which can be identified with the parseOut function.  It returns a list
    #                 containing data for whole genome duplication, fraction genome altered, and fraction LOH,
    #                 or False when the file is unreadable, lacks a needed column, or has no row for that cval.
    #                 Returns: [wgd,fga,frac_loh]
    ######################
    @staticmethod
    def parseSampleQC(qc_file, expected_cval):
        try:
            qc_df = pd.read_csv(qc_file, sep="\t", low_memory=False)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
            return False
        #Sometimes there are missing values, we need to handle those properly.
        #I.E. there is no fga column in some older refit qc.txt files.
        for col in ("fga", "cval", "frac_loh", "wgd"):
            if col not in qc_df.columns:
                return False
        try:
            want = int(float(expected_cval))
        except (TypeError, ValueError):
            return False
        #We want the row with a cval matching our target sample/fit.
        rows = qc_df.loc[pd.to_numeric(qc_df["cval"], errors="coerce") == want]
        if rows.empty:
            return False
        row = rows.iloc[0]
        try:
            return [_to_bool(row["wgd"]), float(row["fga"]), float(row["frac_loh"])]
        except (TypeError, ValueError):
            return False


    ######################
    # parseFacetsQC:  This function accepts the directory map list for a run (a master_file_dict value).
    #                 It returns the sample-level facets_qc pass flag (True/False) for that run's fit, matching
    #                 facets_qc.txt rows on fit_name (which is the nested name for ultra fits, e.g. default/ultra).
    #                 Returns None when the file is unreadable, lacks columns, or has no row for the fit.
    ######################
    @staticmethod
    def parseFacetsQC(sample_dir_map):
        try:
            qc_df = pd.read_csv(sample_dir_map[MetaDictMap.FACETS_QC_FILE], sep="\t", low_memory=False)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError, IndexError, TypeError):
            return None
        if 'fit_name' not in qc_df.columns or 'facets_qc' not in qc_df.columns:
            return None
        names = qc_df['fit_name'].astype(str)
        fit_name = ""
        if len(sample_dir_map) > MetaDictMap.FIT_STATUS:
            parts = str(sample_dir_map[MetaDictMap.FIT_STATUS]).split("#")
            if len(parts) >= 2:
                fit_name = parts[1]
        rows = qc_df.loc[names == fit_name] if fit_name else qc_df.iloc[0:0]
        if rows.empty:
            #Fall back to the fit directory's basename (pre-0.6 behavior).
            base = os.path.basename(os.path.dirname(str(sample_dir_map[MetaDictMap.OUT_FILE])))
            rows = qc_df.loc[names == base]
        if rows.empty:
            return None
        return _to_bool(rows.iloc[0]['facets_qc'])


    ######################
    # parseCNCF:      This function accepts a CNCF file path and parses it for data.
    #                 It will build and return a list of FacetsSegment type objects representing the CNCF.
    ######################
    @staticmethod
    def parseCNCF(cncf_file,adjseg):
        try:
            cncf_df = pd.read_csv(cncf_file, sep="\t", low_memory=False)
            seg_df = pd.read_csv(adjseg, sep="\t", low_memory=False)



            #Extract the relevant columns.
            chrom       = cncf_df['chrom'].tolist()
            loc_start   = cncf_df['loc.start'].tolist()
            loc_end     = cncf_df['loc.end'].tolist()
            cnlr_median = cncf_df['cnlr.median'].tolist()

            cf_em       = cncf_df['cf.em'].tolist()
            cf_base     = cncf_df['cf'].tolist()
            tcn_em      = cncf_df['tcn.em'].tolist()
            lcn_em      = cncf_df['lcn.em'].tolist()
            num_mark    = cncf_df['num.mark'].tolist()
            nhet        = cncf_df['nhet'].tolist()
            mafR        = cncf_df['mafR'].tolist()
            segclust    = cncf_df['segclust'].tolist()
            cnlr_median_clust     = cncf_df['cnlr.median.clust'].tolist()
            mafR_clust  = cncf_df['mafR.clust'].tolist()
            # print(seq_df['seg.mean'].tolist())
            adj_mean    = seg_df['seg.mean'].tolist()


            #Build a FacetsSegment object for each row.
            seg_list = []
            for i in range(len(chrom)):
                cur_seg = FacetsSegment(chrom[i],
                                        loc_start[i],
                                        loc_end[i],
                                        cnlr_median[i],
                                        cf_em[i],
                                        cf_base[i],
                                        tcn_em[i],
                                        lcn_em[i],
                                        num_mark[i],
                                        nhet[i],
                                        mafR[i],
                                        segclust[i],
                                        cnlr_median_clust[i],
                                        mafR_clust[i],
                                        adj_mean[i]
                                        )
                seg_list.append(cur_seg)

            return seg_list

        except (OSError, KeyError, ValueError, TypeError, IndexError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            raise ParseError("Error processing clinical data (parseCNCF). Terminating execution." + ": " + str(e)) from e
    ######################
    # buildFacetsRun:  This function accepts a run key (master_file_dict key) and parses all relevant files and assembles
    #                     all desired values.  It creates and returns an object of type FacetsRun, or False when the
    #                     run's files are empty/inconsistent (counted by the caller as a build failure).
    ######################
    @staticmethod
    def buildFacetsRun(sample_id, facets_metadata):
        cur_dir_map = facets_metadata.master_file_dict.get(sample_id)
        if cur_dir_map is None:
            raise FacetsAPIError("FacetsRun.buildFacetsRun(): %s is not in master_file_dict." % sample_id)

        def warn(message):
            if facets_metadata.run_verbose:
                print(bcolors.WARNING + "\t\tWarning: " + sample_id + " " + message + bcolors.ENDC)

        #If a file exists, but is empty, return a failure.
        for slot, label in ((MetaDictMap.OUT_FILE, "out_file"), (MetaDictMap.CNCF_FILE, "cncf_file"),
                            (MetaDictMap.QC_FILE, "qc_file"), (MetaDictMap.FACETS_QC_FILE, "facets_qc_file"),
                            (MetaDictMap.GENE_FILE, "gene-level file")):
            path = cur_dir_map[slot]
            if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
                warn(label + " is missing or empty!")
                return False

        #Get data from our relevant input files.
        mode           = facets_metadata.hisens_vs_purity
        out_data       = FacetsRun.parseOut(cur_dir_map[MetaDictMap.OUT_FILE], mode)
        if out_data is False:
            warn("out_file could not be read!")
            facets_metadata.failed_samples.append(sample_id)
            return False
        sample_qc_data = FacetsRun.parseSampleQC(cur_dir_map[MetaDictMap.QC_FILE], out_data[1])
        facets_qc_out  = FacetsRun.parseFacetsQC(cur_dir_map)

        #If we were missing information in any of our files, we don't want to proceed with this sample.
        if facets_qc_out is None or facets_qc_out == -1:
            warn("facets_qc_file parse failed when building FacetsRun object! DirMap: " + str(cur_dir_map))
            facets_metadata.failed_samples.append(sample_id)
            return False
        if sample_qc_data is False:
            warn("sample qc parse failed when building FacetsRun object! DirMap: " + str(cur_dir_map))
            facets_metadata.failed_samples.append(sample_id)
            return False

        #Get data from our data mappings. Clinical maps are keyed by the 17-character tumor sample id.
        short_id = short_sample_id(sample_id)
        cur_cancer_type        = facets_metadata.cancer_type_map.get(short_id)
        cur_cancer_type_detail = facets_metadata.cancer_type_detail_map.get(short_id)
        cur_clin_purity        = facets_metadata.clinical_purity_map.get(short_id)
        cur_onko_code          = facets_metadata.onkotree_code_map.get(short_id)
        cur_patient_id         = facets_metadata.patient_id_map.get(short_id)
        cur_tmb                = facets_metadata.cvr_tmb_score_map.get(short_id)
        cur_msi                = facets_metadata.msi_score_map.get(short_id)
        review_status          = (str(cur_dir_map[MetaDictMap.FIT_STATUS]).split("#") + ["", "", ""])[0:3]  #reviewstatus#fit_name#reviewed_by
        cur_review_status, fit_name, reviewed_by = review_status
        fit_class              = cur_dir_map[MetaDictMap.FIT_CLASS] if len(cur_dir_map) > MetaDictMap.FIT_CLASS else "standard"

        #Build the FacetsRun object.
        cur_facets_sample = FacetsRun(sample_id,
                                        cur_patient_id,
                                        cur_dir_map[MetaDictMap.SELECTED_FIT_DIR],
                                        cur_cancer_type,
                                        cur_cancer_type_detail,
                                        out_data[0],
                                        cur_clin_purity,
                                        cur_onko_code,
                                        out_data[2],
                                        out_data[3],
                                        out_data[1],
                                        sample_qc_data[0],
                                        sample_qc_data[1],
                                        sample_qc_data[2],
                                        facets_qc_out,
                                        cur_tmb,
                                        cur_msi,
                                        mode,
                                        cur_review_status,
                                        fit_name,
                                        reviewed_by,
                                        fit_class=fit_class
                                        )

        #Now build in the segments from the CNCF file.
        segments = FacetsRun.parseCNCF(cur_dir_map[MetaDictMap.CNCF_FILE], cur_dir_map[MetaDictMap.ADJUSTED_SEG_FILE])
        for seg in segments:
            if str(seg.cf) == "nan" or str(seg.cf_base) == "nan":
                continue
            cur_facets_sample.addSegment(seg)

        #Now build in the gene-level data.
        geneLevel = FacetsRun.parseGeneLevel(cur_dir_map[MetaDictMap.GENE_FILE])
        for curGene in geneLevel:
            cur_facets_sample.addGene(FacetsGene(*curGene[0:21]))

        return cur_facets_sample

from heapq import merge
from random import sample
import sys
import os.path
import os
import glob
import pandas as pd
import numpy as np
import pickle

import statistics
from datetime import datetime,time
import matplotlib.pyplot as plt



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
# ExtTools:    This class functions as a wrapper to execute certain functions from external sources,
######################
class ExtTools:
    def __init__(self, ref_meta):
        if not isinstance(ref_meta, FacetsMeta):
            print (bcolors.FAIL)
            print ("\t\tError in ExtTools(). ref_meta must be of type FacetsMeta.")
            print (bcolors.ENDC)
            sys.exit()
        self.ref_meta = ref_meta

    #This function will module load rVersionString on juno.
    def loadModule(self, moduleToLoad):
        try:
            os.system("module load " + moduleToLoad)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in ExtTools.load_R_version(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will make a merged file for any file type in the MetaDictMap with a shared header.
    #fileType should be a MetaDictMap.FILETYPE. 
    #This will return a path to the written file.
    def makeMergedFile(self, fileType, outfile_path):
        try:
            #Keep track and don't duplicate any files. 
            #This can happen when there are alt fits included for the same sample.
            #For instance a seg file merge could merge a default and alt fit, but would duplicate
            #the exact values.
            completed_files = [] 
            with open(outfile_path, 'w') as outfile:
                writeHeader = True
                for item in self.ref_meta.master_file_dict:
                    cur_file = self.ref_meta.master_file_dict.get(item)[fileType]
                    base_file = os.path.basename(cur_file)
                    if base_file in completed_files:
                        continue
                    else:
                        completed_files.append(base_file)
                    
                    #Write the header only for the first file.
                    if writeHeader is True:
                        os.system("cat " + cur_file + " > " + outfile_path)
                        writeHeader = False
                    else:
                        os.system("tail -n +2 " + cur_file + " >> " + outfile_path)

            return outfile_path

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in ExtTools.makeMergedFile(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will run ascets.
    #bsub -J "P-0000012-T03-IM3_P-0000012-N01-IM3" -We 1:59 -n 2 bash -c "Rscript run_ascets.R -i /juno/work/ccs/pricea2/data/ascets/input_segs/P-0000012-T03-IM3_P-0000012-N01-IM3_purity.seg -c /juno/work/ccs/pricea2/pipelines/ascets/ascets/genomic_arm_coordinates_hg19.txt -m 0.5 -k F -a 0.7 -o P-0000012-T03-IM3_P-0000012-N01-IM3"
    def runAscets(self, output_dir="", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True):
        try:
            arm_coords_file = ""
            merged_seg_tmp_file = ""

            if output_dir == "":
                print (bcolors.FAIL)
                print ("\t\tError in Tools.runAscets(). No output path provided.")
                print (bcolors.ENDC)
                sys.exit()
            if ref_genome_coords != "hg19" and ref_genome_coords != "hg38":
                print (bcolors.FAIL)
                print ("\t\tError in Tools.runAscets(). ref_genome_coords must be hg19 or hg38.")
                print (bcolors.ENDC)
                sys.exit()
            #Select the appropriate genomic arm coordinate file. 
            else:
                if ref_genome_coords == "hg19":
                    arm_coords_file = os.path.dirname(__file__) + "/../tools/ascets/genomic_arm_coordinates_hg19.txt"
                if ref_genome_coords == "hg38":
                    arm_coords_file = os.path.dirname(__file__) + "/../tools/ascets/genomic_arm_coordinates_hg38.txt"

            #Print some notes on dependencies.
            print("\t\tNote: Ascets requires R 3.6.* and requires tidyverse and data.table libraries.")
            print("\t\tNote: Failure to meet these dependencies will result in an empty ascets run.")

            #Make an output folder if necessary.
            if os.path.isdir(output_dir):
                print ("\t\tAscets folder already exists: " + output_dir)
            else:
                print ("\t\tCreating ascets output folder at: " + output_dir)
                os.system("mkdir " + output_dir)

            #Make an adjusted or unadjusted merged seg file to use as ascets input.
            #This will use everything in the FacetsMeta reference object for ExtTools.
            adj_str = ""
            merged_seg_tmp_file = output_dir + "/merged_seg_data.seg"
            if use_adjusted:
                adj_str = "facetsAdjustedSeg"
                self.makeMergedFile(MetaDictMap.ADJUSTED_SEG_FILE, merged_seg_tmp_file)
            else:
                adj_str = "facetsUnadjustedSeg"
                self.makeMergedFile(MetaDictMap.UNADJUSTED_SEG_FILE, merged_seg_tmp_file)

            noise_str = ""
            keep_noise_short = ""
            if keep_noise:
                keep_noise_short = "T"
                noise_str = "keepNoise"
            else:
                keep_noise_short = "F"
                noise_str = "filterNoise"

            out_log_file = output_dir + "/runLog.out"
            err_log_file = output_dir + "/runLog.err"

            #Make an output prefix:
            out_prefix = "FacetsAscetsRun_" + ref_genome_coords + "_breadth" + str(min_arm_breadth) + "_" + noise_str + "_armAltFracThresh" + str(arm_alt_frac_thresh) + "_" + adj_str + "_" + str(date.today())

            run_ascets_script = os.path.dirname(__file__) + "/../tools/ascets/run_ascets.R"

            #Load required version of R.
            self.loadModule("R/R-3.6.3")

            #Submit ascets run to the queue.
            ascets_cmd = "bsub -J " + out_prefix + " -o " + out_log_file + " -e " + err_log_file + " -We 1:59 -n 2 bash -c \"Rscript " + run_ascets_script + " -i " + merged_seg_tmp_file + " -c " + arm_coords_file + " -m " + str(min_arm_breadth) + " -k " + keep_noise_short + " -a " + str(arm_alt_frac_thresh) + " -o " + out_prefix + " -p " + output_dir + "\""
            print("\t\t\tRunning ascets...")
            os.system(ascets_cmd)
            
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in ExtTools.runAscets(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()    

######################
# FPTools:    This class functions as a wrapper to execute certain functions from Facets Preview.
######################
class FPTools:
    genomic_anno_wrapper_script = "/juno/work/ccs/pricea2/pipelines/facetsAPI/R/fp_generate_genomic_anno_wrapper.R"

    def __init__(self, config_file="/juno/work/ccs/bandlamc/git/ccs-cron/impact_facets/config_facets_preview.json", cbio_maf_file="", cbio_nonsigned_maf_file=""):
        self.config_file             = config_file
        self.cbio_maf_file           = cbio_maf_file
        self.cbio_nonsigned_maf_file = cbio_nonsigned_maf_file

    #This function will module load rVersionString on juno.
    def loadModule(self, moduleToLoad):
        try:
            os.system("module load " + moduleToLoad)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FPTools.load_R_version(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function runs the wrapper script, R/fp_generate_genomic_anno_wrapper.R
    #That script takes these arguments in this order:
    #facets_preview_config_file = args[1]
    #tumor_sample = args[2] -- i.e. P-0000280-T05-IM7
    #sample_id = args[3] -- i.e. P-0071250-T01-IM7_P-0071250-N01-IM7
    #sample_path = args[4]
    #cbio_maf_path = args[5]
    #cbio_maf_nonsignedout_path = args[6]
    def runGenerateGenomicAnnotations(self, short_id, long_id, path):
        try:
            run_cmd = "bsub -W 1:59 -n 1 -R \"rusage[mem=8G]\" Rscript " + self.genomic_anno_wrapper_script + " " + \
                    self.config_file + " " + short_id + " " + long_id + " " + path + " " + \
                    self.cbio_maf_file + " " + self.cbio_nonsigned_maf_file
            print(run_cmd)
            os.system(run_cmd)
            #sys.exit()
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FPTools.runGenerateGenomicAnnotations(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will make backup files of existing manifests for every sample in the provided FacetsMeta object.
    def backupManifests(self, ref_meta):
        try:
            if not isinstance(ref_meta, FacetsMeta):
                print (bcolors.FAIL)
                print ("\t\tError in FPTools.backupManifests(). ref_meta must be of type FacetsMeta.")
                print (bcolors.ENDC)
                sys.exit()
            print("Backing up manifest files.")
            for item in ref_meta.master_file_dict:
                curSample   = ref_meta.master_file_dict.get(item)
                curManifest = curSample[MetaDictMap.MANIFEST_FILE]
                bakManifest = curManifest + "." + str(date.today()) + ".bak"
                backup_cmd  = "cp " + curManifest + " " + bakManifest
                os.system(backup_cmd)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FPTools.backupManifests(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will make backup files of existing facets qc files for every sample in the provided FacetsMeta object.
    def backupFacetsQC(self, ref_meta):
        try:
            if not isinstance(ref_meta, FacetsMeta):
                print (bcolors.FAIL)
                print ("\t\tError in FPTools.backupFacetsQC(). ref_meta must be of type FacetsMeta.")
                print (bcolors.ENDC)
                sys.exit()
            print("Backing up Facets QC files.")
            for item in ref_meta.master_file_dict:
                curSample   = ref_meta.master_file_dict.get(item)
                curQC = curSample[MetaDictMap.FACETS_QC_FILE]
                bakQC = curQC + "." + str(date.today()) + ".bak"
                backup_cmd  = "cp " + curQC + " " + bakQC
                os.system(backup_cmd)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FPTools.backupManifests(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


######################
# MetaDictMap:    This class is a map of files for the master_file_dict in FacetsMeta for better readability. 
#                 When calling master_file_dict, instead of numbered indeces, use MetaDictMap.FILE_NAME.
######################
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
    SAMPLE_BASE_DIR = 12

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

    def __init__(self, clinical_sample_file = "", facets_repo_path = "", hisens_vs_purity="purity", persist_data="no"):
        #Relevant path and file data.
        self.clinical_sample_file    = clinical_sample_file
        self.facets_repo_path        = facets_repo_path
        self.hisens_vs_purity        = hisens_vs_purity
        self.persist_data            = persist_data

        self.hasClinicalMeta = False
        if clinical_sample_file == "":
            self.hasClinicalMeta = False
        else:
            self.hasClinicalMeta = True

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
        self.fit_map                = {} # A map of id -> best/acceptable/default fit.

    def setVerbose(self, doVerbose):
        try:
            run_verbose = doVerbose
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsMeta.setVerbose(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will build the FacetsMeta object. It handles executing parseClinicalSample() and data persistence.
    def buildFacetsMeta(self):
        #Making a meta is going to be the first thing in FACETS API, so this is where we can print the logo.
        self.printLogo()

        #Fill in most of these structures with the parseClinicalSample function.
        #If we are using data persistance, we want to load that from an existing structured file or save the initial parse to an existing file.
        #That way we don't need to process the clinical sample file each run.
        if self.persist_data != "no":
            if os.path.isfile(self.persist_data):
                print("\tLoading previously processed FacetsMeta object.")
                storedMeta = open(self.persist_data, 'rb')
                savedInstance = pickle.load(storedMeta) 
                for k in savedInstance.__dict__.keys():
                    setattr(self, k, getattr(savedInstance, k))
            else:
                self.parseClinicalSample()
                print("\tStoring this FacetsMeta object for future loading.")
                dataToStore = open(self.persist_data, 'wb')
                pickle.dump(self, dataToStore)
        else:
            self.parseClinicalSample()

    #This function will set this FacetsMeta object to select a single run per sample. If allowDefaults is set to True, 
    # default fits will be selected when best/acceptable fits are not available.  If false, only best/acceptable fits will be selected.
    def setSingleRunPerSample(self, doSetSingle, allowDefaults):
        try:
            if not isinstance(doSetSingle, bool):
                print (bcolors.FAIL)
                print ("\t\tError in FacetsMeta.setSingleRunPerSample(). doSetSingle expected True or False value.")
                print (bcolors.ENDC)
                sys.exit()
            if not isinstance(allowDefaults, bool):
                print (bcolors.FAIL)
                print ("\t\tError in FacetsMeta.setSingleRunPerSample(). allowDefaults expected True or False value.")
                print (bcolors.ENDC)
                sys.exit()
            if doSetSingle == False and allowDefaults == True:
                print (bcolors.FAIL)
                print ("\t\tError in FacetsMeta.setSingleRunPerSample(). allowDefaults can only be True with doSetSingle is also True.")
                print (bcolors.ENDC)
                sys.exit()

            FacetsMeta.selectSingleRun         = doSetSingle
            FacetsMeta.use_unreviewed_defaults = allowDefaults

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsMeta.setSingleRunPerSample(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function accepts a file with a single sample ID per line and populates this objects
    #samples_from_file list.  If this list is populated and build_from_file_listing is true,
    #When parseClinicalData runs, it will only include samples listed in the provided file.
    def selectSamplesFromFile(self, infile):
        try:
            self.build_from_file_listing = True
            with open(infile) as fp:
                line = fp.readline()
                while line:
                    #print(line.strip())
                    self.samples_from_file.append(line.strip())
                    line = fp.readline()
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsMeta.selectSamplesFromFile(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #Simple method for printing a FACETS API logo.
    @staticmethod
    def printLogo():
        print("~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~") 
        print("|"+bcolors.OKBLUE+" ______   ______     ______     ______     ______   ______    "+bcolors.OKCYAN+"    ______     ______   __    "+bcolors.ENDC+"|") 
        print("|"+bcolors.OKBLUE+"/\  ___\ /\  __ \   /\  ___\   /\  ___\   /\__  _\ /\  ___\   "+bcolors.OKCYAN+"   /\  __ \   /\  == \ /\ \   "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+"\ \  __\ \ \  __ \  \ \ \____  \ \  __\   \/_/\ \/ \ \___  \  "+bcolors.OKCYAN+"   \ \  __ \  \ \  _-/ \ \ \  "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+" \ \_\    \ \_\ \_\  \ \_____\  \ \_____\    \ \_\  \/\_____\ "+bcolors.OKCYAN+"    \ \_\ \_\  \ \_\    \ \_\ "+bcolors.ENDC+"|")
        print("|"+bcolors.OKBLUE+"  \/_/     \/_/\/_/   \/_____/   \/_____/     \/_/   \/_____/  "+bcolors.OKCYAN+"    \/_/\/_/   \/_/     \/_/ "+bcolors.ENDC+"|")
        print("~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~-===-~~")

    ######################
    # parseClinicalSample:  This function will accept a clinical sample file and
    #                         scan each sample's corresponding facets directory.
    #                         the facets_review.manifest file will be read to identify
    #                         the appropriate fit to use in analysis.  A dictionary mapping
    #                         samples to their appropriate directories is built in master_file_dict.
    #                         This function will also build dictionaries mapping sample id to
    #                         cancer type and to patient id.
    ######################
    def parseClinicalSample(self):
        num_best_fits        = 0
        num_acceptable_fits  = 0
        num_default_fits     = 0
        num_multi_normals    = 0
        num_missing_manifest = 0
        num_missing          = 0

        try:
            print("\tProcessing clinical data...")
            if self.facets_repo_path == "":
                print (bcolors.FAIL)
                print("Error in FacetsMeta.parseClinicalSample(): Facets repo path is not supplied. ")
                print (bcolors.ENDC)
                sys.exit()

            #If there is no clinical sample file, 
            #we just want to build the FacetsMeta object by browsing the directory stucture.
            if self.clinical_sample_file == "":
                startdir=self.facets_repo_path
                print(bcolors.WARNING + "\t\tWarning: No data clinical sample file was provided. ")
                print("\t\t\tUsing FACETS directory structure to build FacetsMeta object: " + startdir + " " + bcolors.ENDC)
                target_ids = []
                for item in os.listdir(startdir):
                    if os.path.isdir(startdir+item):
                        for samplefolder in os.listdir(startdir+item):
                            target_ids.append(samplefolder)
            else:
                clinical_df = pd.read_csv(self.clinical_sample_file, sep="\t", low_memory=False)

                #Extract DMP id's for everything in our clinical sample file.
                print("\t\tIdentifying samples from impact FACETS repository.")

                target_ids           = clinical_df['SAMPLE_ID'].tolist()
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
            if FacetsMeta.selectSingleRun:
                print("\t\tBuilding directory map targeting a best runs for each sample.")
            else:
                print("\t\tBuilding directory map including all runs for each sample.")

            if self.build_from_file_listing == True:
                print("\t\tSelecting samples only from the provided input sample list file.")

            for id in target_ids:
                cur_short_id = id[0:7]
                cur_tumor_id = id[0:17]
                cur_id_dirs = glob.glob(self.facets_repo_path + cur_short_id + "/" + id + "*")

                #print(":::cur_id: " + cur_tumor_id)

                #If we are parsing data from a file, only include the samples we've already read in.
                if self.build_from_file_listing == True:
                    if cur_tumor_id not in self.samples_from_file:
                        continue

                #Keep track of how many samples had the same tumor/id, but multiple normals.
                if len(cur_id_dirs) > 1:
                    num_multi_normals = num_multi_normals + 1

                if not cur_id_dirs:
                    continue
                else:
                    id_with_normal = cur_id_dirs[-1].split('/')[-1]

                    cur_manifest_file  = cur_id_dirs[-1] + "/facets_review.manifest"
                    cur_facets_qc_file = cur_id_dirs[-1] + "/facets_qc.txt" 
                    cur_base_dir       = cur_id_dirs[-1]

                    #If there is no manifest file, skip the sample.
                    if not os.path.exists(cur_manifest_file):
                        num_missing_manifest = num_missing_manifest + 1
                        continue

                    cur_manifest = pd.read_csv(cur_manifest_file, skiprows=1, sep="\t", low_memory=False)
                    manifest_data = cur_manifest[['path', 'review_status', 'fit_name', 'reviewed_by']]

                    #If we are looking for individual fits we will want to look for the best and consider default fits.
                    if FacetsMeta.selectSingleRun:
                        #Grab important data for best/acceptable/default fits.
                        target_best_fit = manifest_data.loc[manifest_data['review_status'].isin(["reviewed_best_fit"])]
                        target_accept_fit = manifest_data.loc[manifest_data['review_status'].isin(["reviewed_acceptable_fit"])]
                        # print(target_best_fit)
                        if(FacetsMeta.use_unreviewed_defaults):
                            target_default_fit = manifest_data.loc[manifest_data['fit_name'].isin(["default"])]
                        # print(target_default_fit)
                        selected_fit_dir = ""
                        #If there is a best fit available, we want to use that.
                        if(not target_best_fit.empty):
                            best_fit_list = target_best_fit.values.tolist()[0]
                            # print(best_fit_list)
                            if best_fit_list[0][-1] != '/':
                                best_fit_list[0] = best_fit_list[0] + "/"
                            selected_fit_dir  = best_fit_list[0] + best_fit_list[2] + "/"

                            curr_fit_status   =  best_fit_list[1] + "#" + best_fit_list[2] + "#" + str(best_fit_list[3])
                            num_best_fits = num_best_fits + 1     
                            self.fit_map[id] = ["Best", selected_fit_dir]    

                        #If there is no best fit, check for acceptable fits.
                        elif(not target_accept_fit.empty):
                            accept_fit_list = target_accept_fit.values.tolist()[0]
                            if accept_fit_list[0][-1] != '/':
                                accept_fit_list[0] = accept_fit_list[0] + "/"
                            selected_fit_dir  = accept_fit_list[0] + accept_fit_list[2] + "/"
                            curr_fit_status   = accept_fit_list[1] + "#" + accept_fit_list[2] + "#" + str(accept_fit_list[3])

                            num_acceptable_fits = num_acceptable_fits + 1
                            self.fit_map[id] = ["Acceptable", selected_fit_dir]
                        #If we have nothing still, and are accepting default fits, use that.
                        elif(FacetsMeta.use_unreviewed_defaults):
                            default_fit_list = target_default_fit.values.tolist()[0]
                            # print(default_fit_list)
                            if default_fit_list[0][-1] != '/':
                                default_fit_list[0] = default_fit_list[0] + "/"
                            selected_fit_dir  = default_fit_list[0] + default_fit_list[2] + "/"
                            num_default_fits = num_default_fits + 1

                            curr_fit_status = default_fit_list[1] + "#" + default_fit_list[2] + "#" + str(default_fit_list[3])
                            self.fit_map[id] = ["Default", selected_fit_dir]

                        #Nothing we want was available, move on from this sample.
                        else:
                            continue

                        #Find the correct file paths depending on if we want hisens or purity.
                        out_file  = ""
                        cncf_file = ""
                        qc_file   = ""
                        cur_out   = ""
                        cur_cncf  = ""
                        gene_level_file = ""
                        adjseg_file = ""
                        if self.hisens_vs_purity == "purity":
                            cur_out      = glob.glob(selected_fit_dir + "/*_purity.out")
                            cur_cncf     = glob.glob(selected_fit_dir + "/*_purity.cncf.txt")
                            cur_adjseq   = glob.glob(selected_fit_dir + "/*_purity_diplogR.adjusted.seg")
                            cur_unadjseq = glob.glob(selected_fit_dir + "/*_purity_diplogR.unadjusted.seg")
                        elif self.hisens_vs_purity == "hisens":
                            cur_out      = glob.glob(selected_fit_dir + "/*_hisens.out")
                            cur_cncf     = glob.glob(selected_fit_dir + "/*_hisens.cncf.txt")
                            cur_adjseq   = glob.glob(selected_fit_dir + "/*_hisens_diplogR.adjusted.seg")
                            cur_unadjseq = glob.glob(selected_fit_dir + "/*_hisens_diplogR.unadjusted.seg")
                        else:
                            print("Error: hisens_vs_purity value should be 'hisens' or 'purity'.")
                            sys.exit()

                        cur_qc         = glob.glob(selected_fit_dir + "*qc.txt")
                        cur_gene_level = glob.glob(selected_fit_dir + "*gene_level.txt")
                        cur_ccf_mafs   = glob.glob(selected_fit_dir + "*ccf.maf")
                        cur_ccf = ""
                        cur_nonsignedout = ""
                        for ccf_maf in cur_ccf_mafs:
                            if 'nonsignedout' in ccf_maf:
                                cur_nonsignedout = ccf_maf
                            else:
                                cur_ccf = ccf_maf
                        
                        #If critical files are missing, skip the sample.
                        if not cur_out or not cur_cncf or not cur_qc or not cur_gene_level or not cur_adjseq or not cur_unadjseq:
                            num_missing = num_missing + 1
                            continue
                        else:
                            out_file        = cur_out[0]
                            cncf_file       = cur_cncf[0]
                            qc_file         = cur_qc[0]
                            gene_level_file = cur_gene_level[0]
                            adjseg_file     = cur_adjseq[0]
                            unadjseg_file   = cur_unadjseq[0]

                        self.long_id_map[id]   = [id_with_normal]

                        self.master_file_dict[id] = [out_file, cncf_file, qc_file, cur_facets_qc_file, selected_fit_dir, gene_level_file, adjseg_file, cur_ccf, cur_nonsignedout, unadjseg_file, cur_manifest_file, curr_fit_status, cur_base_dir]

                        #print(str([out_file, cncf_file, qc_file, cur_facets_qc_file, selected_fit_dir, gene_level_file, adjseg_file]))
                        #break
                    #If we want to read in all fits for each sample, we need to iterate the manifest and build each one out.
                    else:
                        cur_run_list = []
                        for index, row in manifest_data.iterrows():
                            if row['fit_name'] == "Not selected":
                                continue

                            long_id = id_with_normal + "#" + row['fit_name']
                            out_file  = ""
                            cncf_file = ""
                            qc_file   = ""
                            cur_out   = ""
                            cur_cncf  = ""
                            gene_level_file = ""
                            adjseg_file = ""

                            #Handle cases where the path doesn't properly end with a / character.
                            cur_sample_folder = row['path']
                            if cur_sample_folder[-1] != '/':
                                cur_sample_folder = cur_sample_folder + "/"

                            cur_fit_folder = cur_sample_folder + row['fit_name'] + "/"

                            if row['fit_name'] == "reviewed_best_fit":
                                self.fit_map[id] = ["Best", cur_fit_folder]
                            if row['fit_name'] == "reviewed_acceptable_fit":
                                self.fit_map[id] = ["Acceptable", cur_fit_folder]
                            if row['fit_name'] == "default":
                                self.fit_map[id] = ["Default", cur_fit_folder]

                            if self.hisens_vs_purity == "purity":
                                cur_out    = glob.glob(cur_fit_folder + "*_purity.out")
                                cur_cncf   = glob.glob(cur_fit_folder + "*_purity.cncf.txt")
                                cur_adjseq = glob.glob(cur_fit_folder + "/*_purity_diplogR.adjusted.seg")
                                cur_unadjseq = glob.glob(cur_fit_folder + "/*_purity_diplogR.unadjusted.seg")
                            elif self.hisens_vs_purity == "hisens":
                                cur_out    = glob.glob(cur_fit_folder + "*_hisens.out")
                                cur_cncf   = glob.glob(cur_fit_folder + "*_hisens.cncf.txt")
                                cur_adjseq = glob.glob(cur_fit_folder + "/*_hisens_diplogR.adjusted.seg")
                                cur_unadjseq = glob.glob(cur_fit_folder + "/*_hisens_diplogR.unadjusted.seg")
                            else:
                                print("Error: hisens_vs_purity value should be 'hisens' or 'purity'.")
                                sys.exit()

                            cur_qc         = glob.glob(cur_fit_folder + "*qc.txt")
                            cur_gene_level = glob.glob(cur_fit_folder + "*gene_level.txt")
                            cur_ccf_mafs   = glob.glob(cur_fit_folder + "*ccf.maf")
                            cur_ccf = ""
                            cur_nonsignedout = ""
                            for ccf_maf in cur_ccf_mafs:
                                if 'nonsignedout' in ccf_maf:
                                    cur_nonsignedout = ccf_maf
                                else:
                                    cur_ccf = ccf_maf

                            #If critical files are missing, skip the sample.
                            if not cur_out or not cur_cncf or not cur_qc or not cur_gene_level or not cur_adjseq or not cur_unadjseq:
                                num_missing = num_missing + 1
                                continue
                            else:
                                out_file        = cur_out[0]
                                cncf_file       = cur_cncf[0]
                                qc_file         = cur_qc[0]
                                gene_level_file = cur_gene_level[0]
                                adjseg_file     = cur_adjseq[0]
                                unadjseg_file   = cur_unadjseq[0]
                                curr_fit_status   = row["review_status"] + "#" + row["fit_name"] + "#" + str(row["reviewed_by"])

                            cur_run_list.append(long_id)
                            self.master_file_dict[long_id] = [out_file, cncf_file, qc_file, cur_facets_qc_file, cur_fit_folder, gene_level_file, adjseg_file, cur_ccf, cur_nonsignedout, unadjseg_file, cur_manifest_file, curr_fit_status, cur_sample_folder]

                            #print(str([out_file, cncf_file, qc_file, cur_facets_qc_file, cur_fit_folder, gene_level_file, adjseg_file]))
                            #sys.exit()
                        self.long_id_map[id] = cur_run_list

            #print(self.master_file_dict)
            #sys.exit()

            if FacetsMeta.selectSingleRun:
                print("\tSelecting single FacetsRun per FacetsSample.")
                print("\t\tIdentified " + str(num_best_fits) + " best fits.")
                print("\t\tIdentified " + str(num_acceptable_fits) + " acceptable fits.")
                print("\t\tIdentified " + str(num_default_fits) + " default fits.")
            else:
                print("\tSelecting all FacetsRuns for each FacetsSample.")

            print("\tIdentified " + str(num_multi_normals) + " instances of tumor having multiple normals.")
            print("\tSkipped " + str(num_missing) + " samples due to missing files.")
            print("\t\tTotal samples extracted: " + str(len(self.master_file_dict)))

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError processing clinical data (492). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


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
        try:
            print("\t\t\tWriting report to: " + outfile_name)
            with open(outfile_name, 'w') as outfile:
                outfile.write("ID\tfitDir\tfitType\tCancer Type\tCancer Type Detail\tFacets Purity\tClinical Purity\tOnkoCode\tPloidy\tdipLogR\tcVal\tWGD\tFGA\tfrac_loh\tfacets_qc_pass\ttmb\tmsi" + "\n")
                found_list = []
                for cur_run in self.runList:
                    cur_tumor_id = cur_run.id[0:17]
                    cur_fit_item = self.ref_facetsMeta.fit_map.get(cur_tumor_id)

                    outfile.write(str(cur_run.id) + "\t")
                    outfile.write(str(cur_run.fitDir) + "\t")
                    outfile.write(str(cur_fit_item[0]) + "\t")
                    outfile.write(str(cur_run.cancerType) + "\t")
                    outfile.write(str(cur_run.cancerTypeDetail) + "\t")
                    outfile.write(str(cur_run.purity) + "\t")
                    outfile.write(str(cur_run.clinicalPurity) + "\t")
                    outfile.write(str(cur_run.onkoCode) + "\t")
                    outfile.write(str(cur_run.ploidy) + "\t")
                    outfile.write(str(cur_run.dipLogR) + "\t")
                    outfile.write(str(cur_run.cval) + "\t")
                    outfile.write(str(cur_run.wgd) + "\t")
                    outfile.write(str(cur_run.fga) + "\t")
                    outfile.write(str(cur_run.frac_loh) + "\t")
                    outfile.write(str(cur_run.facets_qc) + "\t")
                    outfile.write(str(cur_run.tmb) + "\t")
                    outfile.write(str(cur_run.msi) + "\t")
                    outfile.write("\n")
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.writeReport(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()  

    
    #This function outputs a text file containing information that summarizes a facets dataset
    #This includes metrics like number of best/acceptable fits, number of failed samples that couldn't build, number of cancer types, etc
    #The output contains the date and time.  
    def writeDatasetSummary(self):      
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

            curCancerType = run.cancerType
            if curCancerType not in cancer_type_dict:
                cancer_type_dict[curCancerType] = 0
            cancer_type_dict[curCancerType] += 1
            
            curOnkoCode = run.onkoCode
            if curOnkoCode not in oncocode_dict:
                oncocode_dict[curOnkoCode] = 0
            oncocode_dict[curOnkoCode] += 1

            cur_msi = run.msi
            if cur_msi>=1.5:
                msi_count +=1
            msi_li.append(cur_msi)

            if run.wgd:
                wgd_count+=1
            
            purity_li.append(run.purity)
        
        msi_median = statistics.median(msi_li)
        purity_median = statistics.median(purity_li)
        
        date = datetime.now()
        datestr = str(date.isoformat(sep = "_",timespec='minutes')).replace(":","-")
        outstr = ""
        with open("FacetsDataset_summary_{datest}.tsv".format(datest= datestr),'w') as outf:
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
    def createHistogram(self,variable,inbins=10,interactive=False):

        try:
            #If variable is comprised of many int or float values
            if variable in { 'purity', 'clinicalPurity',  'ploidy', 'dipLogR', 'cval', 'fga', 'frac_loh', 'tmb' , 'msi'}:
                variableli = []
                for run in self.runList:
                    if getattr(run,variable) == 'NA' or type(getattr(run,variable)) == str:
                        pass
                    else:
                        variableli.append(float(getattr(run,variable)))

                # counts, bins = np.histogram(variableli)

                date = datetime.now()
                datestr = str(date.isoformat(sep = "_",timespec='minutes')).replace(":","-")

                n, bins, patches = plt.hist(x=variableli, bins=inbins, color='#0504aa',
                            alpha=0.7, rwidth=0.85, orientation="horizontal")
                plt.ylabel('Value')
                plt.xlabel('Frequency')

                plt.title("Histogram of {var} in constructed FacetsDataset ({datetim})".format(var=variable,datetim = datestr ))
                # plt.xticks(range(int(min(n)),int(max(n))+1))
                plt.yticks(bins)
                if interactive:
                    pass
                else:
                    plt.savefig('./histograms/FacetsDataset_histogram_{var}_{datetim}.png'.format(var=variable,datetim =datestr), bbox_inches="tight")
                    plt.clf()

            #Variable is strings 
            elif variable in {"cancerType", 'oncoCode','wgd', 'review_status'}:
                count_dict = {}
                for run in self.runList:
                    if type(getattr(run,variable)) == bool:

                        if str(getattr(run,variable)) not in count_dict:
                            count_dict[str(getattr(run,variable))] = 0
                        count_dict[str(getattr(run,variable))] += 1

                    else:
                        if getattr(run,variable) not in count_dict:
                            count_dict[getattr(run,variable)] = 0
                        count_dict[getattr(run,variable)] += 1

                #Make histogram 
                date = datetime.now()
                datestr = str(date.isoformat(sep = "_",timespec='minutes')).replace(":","-")

                plt.barh(list(count_dict.keys()), list(count_dict.values()), color='g')
                plt.xlabel('Value')
                plt.ylabel('Frequency')
                plt.title("Histogram of {var} in constructed FacetsDataset ({datetim})".format(var=variable,datetim = datestr ))
                plt.autoscale()
                if interactive:
                    pass
                else:
                    plt.savefig('./FacetsDataset_histogram_{var}_{datetim}.png'.format(var=variable,datetim =datestr ),bbox_inches="tight")
                    plt.clf()
                
            else:
                print (bcolors.FAIL)
                print ("\t\tError in FacetsDataset.createHistogram(). Enter an appropriate object from FacetsRun")
                print (e)
                print (bcolors.ENDC)

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.createHistogram(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)


    ######################
    # copyDatasetToFolder:  This function will create a copy of a the folders in the facets dataset object
    #                       to a specified output location. Directory structure will be preserved, /all/P-12345/P-12345_T01...
    #                       Note that this function will not overwrite existing files.
    ######################     
    def copyDatasetToFolder(self, output_folder):
        try:
            if output_folder == self.ref_facetsMeta.facets_repo_path:
                print (bcolors.FAIL)
                print ("\t\tError in FacetsDataset.copyDatasetToFolder(). Destination path cannot be the same as the source path.")
                print (e)
                print (bcolors.ENDC)
                sys.exit()
            
            target_base_dir = output_folder + "/all/"
            #If the target directory doesnt exist, create it.
            if not os.path.exists(output_folder):
                os.system("mkdir " + output_folder)
                os.system("chmod -R 775 " + output_folder)
            #Make an /all/ directory to maintain structure.
            if not os.path.exists(target_base_dir):
                os.system("mkdir " + target_base_dir)
                os.system("chmod -R 775 " + target_base_dir)

            #Copy everything to the appropriate location.
            for key in self.sampleList:
                cur_sample = self.sampleList.get(key)
                cur_run = cur_sample.runs[0] #A sample might have multiple runs, but we don't care here, just want the base directory.
                run_base_dir = os.path.dirname(os.path.dirname(cur_run.fitDir)) #Get the long folder path (P-1234542-IM3-T01_P1234542-IM3-N01/)
                run_short_parent_dir = os.path.dirname(run_base_dir) #Get the short folder (P-01234/)
                short_dir_str = run_short_parent_dir.split("/")[-1]
                long_dir_str = run_base_dir.split("/")[-1]
                new_short_parent_dir = target_base_dir + short_dir_str
                if not os.path.exists(new_short_parent_dir):
                    os.system("mkdir " + new_short_parent_dir)
                    os.system("chmod -R 775 " + new_short_parent_dir)
                target_copy_path = new_short_parent_dir + "/" + long_dir_str
                if not os.path.exists(target_copy_path):
                    os.system("cp -R -n " + run_base_dir + " " + target_copy_path)

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.copyDatasetToFolder(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)

        
    #This function will calculate and write calculated purity data for each sample.  If use_base_cf is true, cf will be used
    #rather than cf.em.  
    def writePurityCFs(self, outfile_path, use_base_cf=False, use_long_id=True):
        try:
            with open(outfile_path, 'w') as outfile:
                headerString = "ID\tPurity\n"
                outfile.write(headerString)
                for curSample in self.runList:
                    curPurity = curSample.calculatePurityByCF(use_base_cf)
                    if use_long_id:
                        outfile.write(self.ref_facetsMeta.long_id_map.get(curSample.id)[0] + "\t" + str(curPurity) + "\n")
                    else:
                        outfile.write(curSample.id + "\t" + str(curPurity) + "\n")

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.writePurityCFs(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will run alteration analysis on the FacetsDataset.
    def runAlterationAnalysis(self):
        try:
            for sample in self.sampleList:
                cur_sample = self.sampleList.get(sample)
                for cur_run in cur_sample.runs:
                    cur_run.defineArms()
                    cur_run.defineAllLohAndGains()
                    cur_run.defineCFLevels()
                    cur_run.defineSampleLevelAlterationCall(FacetsRun.alteration_cov_min, FacetsRun.alteration_arm_num)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.runAlterationAnalysis(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will write an output file containing arm level data and alteration calls.
    #Output file will be id, cancer_type, cancer_type_detail, purity, onkoCode, tmb, msi, wgd, fga, frac_loh, isHypoploidy, isGainHeavy
    #                        1p_alteration, 1q_alteration...22q_alteration, 
    #                        1p_percentAltered, 1q_percentAltered...22q_percentAltered
    #                        1p_cf_value, 1q_cf_value...22q_cf_value, 
    #                        1p_cf_armCov, 1q_cf_armCov...22q_cf_value,
    #                        1p_level, 1q_level...22q_level
    def writeAlterationData(self, outfile_path):
        try:
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
                    #print("hi" + curSample.cancerType)
                    if curSample.id == "P-0054843-T01-IM6" or curSample.id == "P-0058941-T01-IM7":
                        #print("SKIP")
                        continue

                    #curSample.printSample()
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
        
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.writeAlterationData(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will print a specific sample in this facets dataset.
    def printFacetsSampleById(self, id):
        if id not in self.sampleList:
            print (bcolors.WARNING)
            print ("\t\Warning in FacetsDataset.printFacetsSampleById(). ID not found: " + str(id))
            print (e)
            print (bcolors.ENDC)
        else:
            target_sample = self.sampleList.get(id)
            target_sample.printSample()


    #This function will build a facets data set using whatever filters are set.
    def buildFacetsDataset(self):
        try:
            if not isinstance(self.ref_facetsMeta, FacetsMeta):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.buildFacetsDataset(): buildFacetsDataset must be an initialized object of type FacetsMeta.")
                print (bcolors.ENDC)
                sys.exit()

            total_samples_prepped = 0
            num_build_failed      = 0
            num_file_failed       = 0
            
            print("\tApplying filtering and building FacetsDataset beginning with " + str(len(self.ref_facetsMeta.master_file_dict)) + " samples.")

            #Go through our master file dictionary and build facets sample objects.
            for key in self.ref_facetsMeta.master_file_dict.keys():
                #print(key)
                if total_samples_prepped % 1000 == 0 and total_samples_prepped != 0:
                    print("\t\tSamples Processed: " + str(total_samples_prepped))

                cur_dir_map = self.ref_facetsMeta.master_file_dict.get(key)

                #CancerType Filter.
                if self.filterCancerType:
                    curCancerType = self.ref_facetsMeta.cancer_type_map.get(key)
                    if curCancerType not in self.selectedCancerTypes:
                        continue

                #CancerTypeDetail Filter.
                if self.filterCancerTypeDetail:
                    curCancerDetailType = self.ref_facetsMeta.cancer_type_detail_map.get(key)
                    if curCancerDetailType not in self.curCancerDetailType:
                        continue

                #Clinical Purity Filter.
                if self.filterClinicalPurity:
                    curClinPurity = self.ref_facetsMeta.clinical_purity_map.get(key)
                    if curClinPurity == "NA":
                        continue
                    curClinPurity = float(curClinPurity)
                    if curClinPurity < self.minClinPurity or curClinPurity > self.maxClinPurity:
                        continue

                #OnkoCode Filter.
                if self.filterOnkoCode:
                    curOnkoCode = self.ref_facetsMeta.onkotree_code_map.get(key)
                    if curOnkoCode not in self.selectedOnkoCodes:
                        continue

                #TMB Filter.
                if self.filterTMB:
                    cur_tmb = self.ref_facetsMeta.cvr_tmb_score_map.get(key)
                    if cur_tmb == "NA":
                        continue
                    cur_tmb = float(cur_tmb)
                    if cur_tmb < self.minTMB or cur_tmb > self.maxTMB:
                        continue

                #MSI Filter.
                if self.filterMSI:
                    cur_msi = self.ref_facetsMeta.msi_score_map.get(key)
                    if cur_msi == "NA":
                        continue
                    cur_msi = float(cur_msi)
                    if cur_msi < self.minMSI or cur_msi > self.maxMSI:
                        continue

                #Only parse the .out file once if necessary, not for each individual filter.
                cur_purity  = 0
                cur_cVal    = 0
                cur_ploidy  = 0
                cur_dipLogR = 0
                out_data    = []
                
                if self.filterPurity or self.filterPloidy or self.filterDipLogR or self.filterCval or self.filterWGD or self.filterWGD or self.filterFGA or self.filterFracLoH:
                    out_data    = FacetsRun.parseOut(cur_dir_map[0])
                    if out_data is False:
                        num_file_failed += 1
                        continue
                    cur_purity  = out_data[0]
                    cur_cVal    = out_data[1]
                    cur_ploidy  = out_data[2]
                    cur_dipLogR = out_data[3]

                #Only parse the sampleQC file once if necessary, not for each individual filter.
                cur_wgd      = 0
                cur_fga      = 0
                cur_frac_loh = 0
                if self.filterWGD or self.filterFGA or self.filterFracLoH:
                    sample_qc_data = FacetsRun.parseSampleQC(cur_dir_map[2], out_data[1])
                    if sample_qc_data is False:
                        num_file_failed += 1
                        continue
                    cur_wgd      = sample_qc_data[0]
                    cur_fga      = sample_qc_data[1]
                    cur_frac_loh = sample_qc_data[2]

                #Purity Filter.
                if self.filterPurity:
                    if cur_purity == "NA":
                        continue
                    cur_purity = float(cur_purity)
                    #Note: This is some strange nuance of Facets that these specific values can't be trusted.  Ask Allison.
                    if cur_purity == 0.3 or cur_purity == 0.30:
                        continue
                    if cur_purity < self.minPurity or cur_purity > self.maxPurity:
                        continue

                #Ploidy Filter.
                if self.filterPloidy:
                    if cur_ploidy == "NA":
                        continue
                    cur_ploidy = float(cur_ploidy)
                    if cur_ploidy < self.minPloidy or cur_ploidy > self.maxPloidy:
                        continue

                #DipLogR Filter.
                if self.filterDipLogR:
                    if cur_dipLogR == "NA":
                        continue
                    cur_dipLogR = float(cur_dipLogR)
                    if cur_dipLogR < self.minDipLogR or cur_dipLogR > self.maxDipLogR:
                        continue

                #cVal Filter.
                if self.filterCval:
                    if cur_cVal == "NA":
                        continue
                    cur_cVal = float(cur_cVal)
                    if cur_cVal < self.minC or cur_cVal > self.maxC:
                        continue

                #WGD Filter.
                if self.filterWGD:
                    #If we want to keep WGD samples, we want to ignore when cur_wgd is false.
                    if self.keepWGD:
                        if not cur_wgd:
                            continue
                    #If we want to get rid of WGD samples, we want to ignore when cur_wgd is true.
                    else:
                        if cur_wgd:
                            continue

                #FGA Filter.
                if self.filterFGA:
                  if cur_fga < self.minFGA or cur_fga > self.maxFGA:
                        continue

                #Frac LoH Filter.
                if self.filterFracLoH:
                    if cur_frac_loh < self.minLoh or cur_frac_loh > self.maxLoh:
                        continue

                #Facets QC Filter.
                if self.filterFacetsQC:
                    facets_qc_out = FacetsRun.parseFacetsQC(cur_dir_map)
                    if facets_qc_out == -1:
                        num_file_failed += 1
                        continue
                    #If we want to keep facets QC Pass samples, we want to ignore when facets_qc_out is false.
                    if self.keepQCPass:
                        if not facets_qc_out:
                            continue
                    #If we want to get rid of facets QC Pass samples, we want to ignore when facets_qc_out is true.
                    else:
                        if facets_qc_out:
                            continue

                #We've passed all filters, we can build the FacetsRun and FacetsSample for this sample.
                cur_run = FacetsRun.buildFacetsRun(key,self.ref_facetsMeta)
                if cur_run is False:
                    num_build_failed += 1
                    continue
                else:
                    cur_sample_id = ""
                    if FacetsMeta.selectSingleRun:
                        #print("SingleRun - " + str(key))
                        cur_sample_id = key
                        if cur_sample_id in self.ref_facetsMeta.long_id_map:
                            cur_long_id = self.ref_facetsMeta.long_id_map.get(cur_sample_id)[0]
                            if cur_long_id in self.sampleList:
                                cur_sample = self.sampleList.get(cur_long_id)
                                cur_sample.addRun(cur_run)
                            else:
                                cur_sample = FacetsSample(cur_long_id)
                                cur_sample.addRun(cur_run)
                                self.sampleList[cur_long_id] = cur_sample

                                #cur_sample.printSample()
                                #sys.exit()
                        else:
                            print (bcolors.FAIL)
                            print ("\t\tError in FacetsDataset.buildFacetsDataset(). No long id found for " + str(key))
                            print (e)
                            print (bcolors.ENDC)
                            sys.exit()
                    else:
                        #print("Multirun - " + str(key))
                        #When we use all runs in a sample, we are using long id (ID_runID, i.e. P-0000067-T01-IM3_P-0000067-N01-IM3#altDiplogR-02).
                        #We are splitting for now on that # character to get the short id.
                        cur_sample_id = key.split('#')[0]
                        if cur_sample_id in self.sampleList:
                            cur_sample = self.sampleList.get(cur_sample_id)
                            cur_sample.addRun(cur_run)
                        else:
                            cur_sample = FacetsSample(cur_sample_id)
                            cur_sample.addRun(cur_run)
                            self.sampleList[cur_sample_id] = cur_sample
                    self.runList.append(cur_run)
                    total_samples_prepped = total_samples_prepped + 1

            self.numSamples = len(self.runList)
            self.num_file_failed  = num_file_failed
            self.num_build_failed = num_build_failed
            print("\tFacetsDataset built with " + str(self.numSamples) + " samples.")
            print("\t\t" + str(num_file_failed) + " samples failed due to file formatting errors.")
            print("\t\t" + str(num_build_failed) + " samples failed during FacetsRun object construction.")

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.buildFacetsDataset(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function sets a cancer type filter.  The selectedCancerTypes parameter should be a list of strings.
    def setCancerTypeFilter(self,selectedCancerTypes):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.selectedCancerTypes(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not isinstance(selectedCancerTypes, list):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.selectedCancerTypes(): selectedCancerTypes must be a list of string values. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterCancerType    = True
                self.selectedCancerTypes = selectedCancerTypes
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.selectedCancerTypes(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a cancer type detail filter.  The selectedCancerDetailTypes parameter should be a list of strings.
    def setCancerTypeDetailFilter(self,selectedCancerDetailTypes):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCancerTypeDetailFilter(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not isinstance(selectedCancerDetailTypes, list):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCancerTypeDetailFilter(): selectedCancerDetailTypes must be a list of string values. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterCancerTypeDetail    = True
                self.selectedCancerDetailTypes = selectedCancerDetailTypes
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setCancerTypeDetailFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a facets purity filter.  
    def setPurityFilter(self, minPurity, maxPurity=1.0):
        try:
            if not (isinstance(minPurity, float) or isinstance(minPurity, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCancerTypeDetailFilter(): minPurity must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxPurity, float) or isinstance(maxPurity, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCancerTypeDetailFilter(): maxPurity must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterPurity = True
                self.minPurity    = minPurity
                self.maxPurity    = maxPurity
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setCancerTypeDetailFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a clinical purity filter.  The selectedCancerDetailTypes parameter should be a list of strings.
    def setClinicalPurityFilter(self, minClinPurity, maxClinPurity=1.0):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setClinicalPurityFilter(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not (isinstance(minClinPurity, float) or isinstance(minClinPurity, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setClinicalPurityFilter(): minClinPurity must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxClinPurity, float) or isinstance(maxClinPurity, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setClinicalPurityFilter(): maxPurity must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterClinicalPurity = True
                self.minClinPurity        = minClinPurity
                self.maxClinPurity        = maxClinPurity
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setClinicalPurityFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a onkoCode filter.  The selectedOnkoCodes parameter should be a list of strings.
    def setOnkoCodeFilter(self, selectedOnkoCodes):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setOnkoCodeFilter(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not isinstance(selectedOnkoCodes, list):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setOnkoCodeFilter(): selectedOnkoCodes must be a list of string values. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterOnkoCode    = True
                self.selectedOnkoCodes = selectedOnkoCodes
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setOnkoCodeFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a ploidy filter.  The min and max ploidy parameters should be a numeric.
    def setPloidyFilter(self, minPloidy, maxPloidy):
        try:
            if not (isinstance(minPloidy, float) or isinstance(minPloidy, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setPloidyFilter(): minPloidy must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxPloidy, float) or isinstance(maxPloidy, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setPloidyFilter(): maxPloidy must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterPloidy = True
                self.minPloidy    = minPloidy
                self.maxPloidy    = maxPloidy
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setPloidyFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a diplogR filter.  The min and max parameters should be a numeric.
    def setDipLogRFilter(self, minDipLogR, maxDipLogR):
        try:
            if not (isinstance(minDipLogR, float) or isinstance(minDipLogR, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setDipLogRFilter(): minDipLogR must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxDipLogR, float) or isinstance(maxDipLogR, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setDipLogRFilter(): maxDipLogR must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterDipLogR = True
                self.minDipLogR    = minDipLogR
                self.maxDipLogR    = maxDipLogR
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setDipLogRFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a cVal filter.  The min and max parameters should be a numeric.
    def setCvalFilter(self, minC, maxC):
        try:
            if not (isinstance(minC, float) or isinstance(minC, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCvalFilter(): minC must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxC, float) or isinstance(maxC, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setCvalFilter(): maxC must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterCval = True
                self.minC       = minC
                self.maxC       = maxC
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setCvalFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a onkoCode filter.  The keepWGD parameter should be a boolean
    def setWGDFilter(self, keepWGD):
        try:
            if not isinstance(keepWGD, bool):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setWGDFilter(): keepWGD must be True or False. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterWGD = True
                self.keepWGD   = keepWGD
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setWGDFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a FGA filter.  The min and max parameters should be a numeric.
    def setFGAFilter(self, minFGA, maxFGA=1.0):
        try:
            if not (isinstance(minFGA, float) or isinstance(minFGA, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setFGAFilter(): minFGA must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxFGA, float) or isinstance(maxFGA, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setFGAFilter(): maxFGA must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterFGA = True
                self.minFGA    = minFGA
                self.maxFGA    = maxFGA
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setFGAFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a FGA filter.  The min and max parameters should be a numeric.
    def setFracLohFilter(self, minLoh, maxLoh=1.0):
        try:
            if not (isinstance(minLoh, float) or isinstance(minLoh, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setFracLohFilter(): minLoh must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxLoh, float) or isinstance(maxLoh, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setFracLohFilter(): maxLoh must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterFracLoH = True
                self.minLoh        = minLoh
                self.maxLoh        = maxLoh
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setFracLohFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a facetsQC filter.  The keepQCPass parameter should be a boolean
    def setFacetsQCFilter(self, keepQCPass):
        try:
            if not isinstance(keepQCPass, bool):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setFacetsQCFilter(): keepQCPass must be True or False. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterFacetsQC = True
                self.keepQCPass     = keepQCPass
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setFacetsQCFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a TMB filter.  The min and max parameters should be numeric.
    def setTMBFilter(self, minTMB, maxTMB):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setTMBFilter(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not (isinstance(minTMB, float) or isinstance(minTMB, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setTMBFilter(): minTMB must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxTMB, float) or isinstance(maxTMB, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setTMBFilter(): maxTMB must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterTMB = True
                self.minTMB    = minTMB
                self.maxTMB    = maxTMB
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setTMBFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function sets a MSI filter.  The min and max parameters should be numeric.
    def setMSIFilter(self, minMSI, maxMSI):
        try:
            if not self.ref_facetsMeta.hasClinicalMeta:
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setMSIFilter(): Clinical data was not provided for this run. ")
                print("Cannot apply this filter without clinical metadata.  Provide a data_clinical_sample file or remove this filter.")
                print (bcolors.ENDC)
                sys.exit()
            if not (isinstance(minMSI, float) or isinstance(minMSI, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setMSIFilter(): minMSI must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            elif not (isinstance(maxMSI, float) or isinstance(maxMSI, int)):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.setMSIFilter(): maxMSI must be a numeric value. ")
                print (bcolors.ENDC)
                sys.exit()
            else:
                self.filterMSI = True
                self.minMSI    = minMSI
                self.maxMSI    = maxMSI
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsDataset.setMSIFilter(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


######################
# FacetsGene:    This class represents a single facets gene as represented in the gene_level.txt file.
######################
class FacetsGene:
    def __init__(self, gene, gene_start, gene_end, seg_start, seg_end, seg_length, cf, tcn, lcn, cn_state, filter, tsg, seg, median_cnlr_seg, segclust, mcn, genes_on_seg, gene_snps, gene_het_snps, spans_segs):
        self.gene            = gene
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
        try:
            if not isinstance(geneToCompare, FacetsGene):
                print (bcolors.FAIL)
                print("Error: FacetsGene.compareSegments(): Can only compare objects of type FacetsGene.")
                print (bcolors.ENDC)
            else:
                if self.gene == geneToCompare.gene and self.gene_start == geneToCompare.gene_start and self.gene_end == geneToCompare.gene_end:
                    return True
                else:
                    return False
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsSegment.compareSegments(): Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


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
        try:
            if not isinstance(segToCompare, FacetsSegment):
                print (bcolors.FAIL)
                print("Error: Can only compare objects of type FacetsSegment.")
                print (bcolors.ENDC)
            else:
                if self.chrom == segToCompare.chrom and self.start == segToCompare.start and self.end == segToCompare.end:
                    return True
                else:
                    return False
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsSegment.compareSegments(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

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
            print (bcolors.FAIL)
            print("Error: FacetsSample.addRun(): Can only add objects of type FacetsRun.")
            print (bcolors.ENDC)
            sys.exit()
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
                 review_status, fit_name, review_author ):
        self.segments         = []
        self.genes            = []
        self.id               = id
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
        


    #Add a segment to the segment array of type FacetsSegment.
    def addSegment(self, seg):
        if not isinstance(seg, FacetsSegment):
            print (bcolors.FAIL)
            print("Error: FacetsRun.addSegment(): Can only add objects of type FacetsSegment.")
            print (bcolors.ENDC)
            sys.exit()
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
                    print (bcolors.FAIL)
                    print("Error in FacetsRun.removeSegment(): Remove failed to shrink segment array.")
                    print (bcolors.ENDC)
                    sys.exit()
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
            print (bcolors.FAIL)
            print("Error: FacetsRun.addGene(): Can only add objects of type FacetsGene.")
            print (bcolors.ENDC)
            sys.exit()
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
                    print (bcolors.FAIL)
                    print("Error in FacetsRun.removeGene(): Remove failed to shrink gene array.")
                    print (bcolors.ENDC)
                    sys.exit()
        return False

    #This function, given a target chromosome and target arm, will return a list of all segments for this sample that fall in the requested regions.
    def getSegmentsByChromAndArm(self, target_chrom, target_arm):
        try:
            if(int(target_chrom) <= 0 or int(target_chrom) >= 23):
                print (bcolors.FAIL)
                print("Error in FacetsRun.getSegmentsByChromAndArm(): target_chrom must be between 1 and 23.")
                print (bcolors.ENDC)
                sys.exit()
            if(target_arm != "p" and target_arm != "q"):
                print (bcolors.FAIL)
                print("Error in FacetsRun.getSegmentsByChromAndArm(): target_arm must be 'p' or 'q'.")
                print (bcolors.ENDC)
                sys.exit()

            arm_to_fetch = str(target_chrom) + str(target_arm)
        
            target_seg_list = []
            for i in range(len(self.segments)):
                #self.segments[i].printSegment()
                if self.segments[i].arm == arm_to_fetch:
                    target_seg_list.append(self.segments[i])

            return target_seg_list 
            
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getSegmentsByChromAndArm(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function, given a target chromosome, will return a list of all segments for this sample that fall in the requested regions.
    def getSegmentsByChrom(self, target_chrom):
        try:
            if(target_chrom <= 0 or target_chrom >= 23):
                print (bcolors.FAIL)
                print("Error in FacetsRun.getSegmentsByChrom(): target_chrom must be between 1 and 23.")
                print (bcolors.ENDC)
                sys.exit()

            target_seg_list = []
            for i in range(len(self.segments)):
                if self.segments[i].chrom == target_chrom:
                    target_seg_list.append(self.segments[i])

            return target_seg_list 
            
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getSegmentsByChrom(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will calculate purity based on all of the segments associated with this sample.
    #Purity is the maximum cf value that is not 1 or NA. if use_base is set to true, this function
    #will calculate the purity based on the non .em cf values.  Otherwise cf.em will be used.
    def calculatePurityByCF(self, use_base=False):
        try:
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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.calculatePurityByCF(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will return the number of arms with LOH over 'percent_lost' in this sample.
    def getNumberLossArms(self, percent_lost):
        try:
            total_lost = 0
            for arm in self.alterations:
                #Alterations: {chr_arm -> [alteration, percentAltered]}
                cur_arm = self.alterations.get(arm)
                if (cur_arm[0] == "LOH" or cur_arm[0] == "LOH/GAIN") and cur_arm[1] >= percent_lost:
                    total_lost = total_lost + 1

            return total_lost

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getNumberLossArms(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will return the number of arms with GAIN over 'percent_gain' in this sample.
    def getNumberGainArms(self, percent_gain):
        try:
            total_gain = 0
            for arm in self.alterations:
                #Alterations: {chr_arm -> [alteration, percentAltered]}
                cur_arm = self.alterations.get(arm)
                if (cur_arm[0] == "GAIN" or cur_arm[0] == "LOH/GAIN") and cur_arm[1] >= percent_gain:
                    total_gain = total_gain + 1

            return total_gain

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getNumberGainArms(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will sort the self.cfLevel_arms into and return an indexed list of arm levels.
    #The function returns two lists, 1 indexed (ignore 0 position), for [[p_arms],[q_arms]].
    def getSortedCFs(self):
        try:
            if len(self.cfLevel_arms) == 0:
                print (bcolors.FAIL)
                print ("\t\tError in FacetsRun.getSortedCFs(). cfLevels is empty. Try running defineCFLevels().")
                print (e)
                print (bcolors.ENDC)
                sys.exit()

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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getSortedCFs(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function returns cf data for a specified arm, i.e. 1p.
    #It is used because arm_cf_data isn't sorted since its a dictionary and sometimes we want to work
    #with data in a specific order, i.e. when writing files.
    def getCFbyArm(self, target_arm):
        try:
            return self.arm_cf_data.get(target_arm)
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getCFbyChromAndArm(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    # This function will calculate a CF value covering the largest segment of a target_arm.
    # I.E. if an arm is 70% cf=.52 and 30% cf=.23, this function will return .52.  
    # match_cf_range is the maximum distance that two segments can be from one another in order to be considered the same CF.
    # This function returns [longest_cf_value, percent_of_arm_covered].
    def getArmLevelCF(self, target_arm, match_cf_range=0.01):
        try: 
            if target_arm == "undefined":
                return False

            curArmSegs = self.getSegmentsByChromAndArm(target_arm[:-1], target_arm[-1])
            #If there are no segments for this arm, we want to return some negative value instead of just crashing.
            if len(curArmSegs) == 0:
                return False

            cf_percents = {} # A map of cf values and how much they cover. i.e: .37 -> .25, .66 -> .65
            for curSeg in curArmSegs:
                if curSeg.arm == "undefined" or curSeg.percentArm == -1:
                    print (bcolors.FAIL)
                    print("Error in FacetsRun.getArmLevelCF(): Undefined arms in " + self.id + ". Run defineArms().")
                    print (bcolors.ENDC)
                    sys.exit()

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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.getArmLevelCF(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

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
        try:
            if len(self.segments) == 0:
                print (bcolors.FAIL)
                print("Error in FacetsRun.defineCFLevels(): No segments available for sample " + self.id)
                print (bcolors.ENDC)
                sys.exit()
            if len(self.alterations) == 0:
                print (bcolors.FAIL)
                print("Error in FacetsRun.defineCFLevels(): Alterations have not been defined for " + self.id + ". Run defineArms() first.")
                print (bcolors.ENDC)
                sys.exit()
            
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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.defineCFLevels(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #Check segments to determine p or q arm. If centromere is crossed, split the segment into two new segments.
    def defineArms(self):
        try:
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
                        print (bcolors.FAIL)
                        print("Error in FacetsRun.defineArms(): End position occurs before start position.")
                        self.segments[i].printSegment()
                        print (bcolors.ENDC)
                        sys.exit()

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
                        print (bcolors.FAIL)
                        print ("\t\tError in FacetsRun.defineArms(). Failed to remove centromere split segment.")
                        print (bcolors.ENDC)
                        sys.exit()

                #Add segments created from splitting old segments.
                for i in range(len(segsToAdd)):
                    addOk = self.addSegment(segsToAdd[i])
                    if not addOk:
                        print (bcolors.FAIL)
                        print ("\t\tError in FacetsRun.defineArms(). Failed to add centromere split segment.")
                        print (bcolors.ENDC)
                        sys.exit()
        
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.defineArms(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This function will check this sample and determine if it is a hypoploidy sample.
    #It will make these determinations based on a minimum percent of arm covered for min_arm_num of samples.
    #This function will set isHypoploidy and isGainHeavy for this sample.
    def defineSampleLevelAlterationCall(self, min_percent_arm_cov, min_arm_num):
        try:
            if self.getNumberLossArms(min_percent_arm_cov) >= min_arm_num:
                self.isHypoploidy = True
            if self.getNumberGainArms(min_percent_arm_cov) >= min_arm_num:
                self.isGainHeavy = True
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.defineHypoploidy(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    #This is a helper function for defineAllLohAndGains().  
    #It will perform LoH and Gain calculations on a set of segments for a given chromosome and arm.
    def calculateLohAndGains(self, selected_segs, target_chr, target_arm):
        try:
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
                print (bcolors.FAIL)
                print ("\t\tError in FacetsRun.calculateLohAndGains(). Total LCN=0 length or TCN>=Gain is greater than arm length.")
                print (bcolors.ENDC)
                sys.exit()
            else:
                total_lcn0_percent    = float((total_lcn0_len / cur_arm_len) * 100.0)
                total_tcnGain_percent = float((total_tcnGain_len / cur_arm_len) * 100.0)

                #Each chr/arm combination must be unique in the alterations map.
                if cur_chr_arm in self.alterations:
                    print (bcolors.FAIL)
                    print ("\t\tError in FacetsRun.calculateLohAndGains(). Duplicate key in alteration map.")
                    print (bcolors.ENDC)
                    sys.exit()

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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.calculateLohAndGains(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()
    

    ######################
    # defineAllLohAndGains:  This function will process all segments for this FacetsRun object and populate
    #                        the 'alterations' map for the sample.  A sample is determined to be LoH if it is 0 LCN for
    #                        more than the defined FacetsSegment threshold percent of of the arm for all segments on that arm.
    #                        A gain is determined by a correspondingly TCN greater than FacetsSegment.min_gain_tcn for FacetsSegment.percent_arm_gain.
    ######################
    def defineAllLohAndGains(self):
        try:
            for i in range(1,23):
                #Handle p and q arms for each chromosome.
                cur_segs_p    = self.getSegmentsByChromAndArm(i, "p")
                cur_segs_q    = self.getSegmentsByChromAndArm(i, "q")
                self.calculateLohAndGains(cur_segs_p, i, "p")
                self.calculateLohAndGains(cur_segs_q, i, "q")

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.defineAllLohAndGains(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


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
    #            It will extract and return purity, cval, ploidy, and dipLogR values as a list.
    #            Returns: [purity,cval,ploidy,dipLogR]
    ######################
    @staticmethod
    def parseOut(out_file):
        try:
            cur_purity  = ""
            cur_cval    = ""
            cur_ploidy  = ""
            cur_dipLogR = ""

            with open(out_file, 'r') as input:
                for line in input:
                    #Get purity value.
                    if '# Purity' in line:
                        if not isinstance(line.split('=')[1].strip(), float):
                            cur_purity = line.split('=')[1].strip()
                        else:
                            cur_purity = float(line.split('=')[1].strip())

                    #Get cval depending on hisens or purity.
                    if FacetsMeta.hisens_vs_purity == "hisens":
                        if '# cval' in line:
                            if not isinstance(line.split('=')[1].strip(), float):
                                cur_cval = line.split('=')[1].strip()
                            else:
                                cur_cval = float(line.split('=')[1].strip())
                    elif FacetsMeta.hisens_vs_purity == "purity":
                        if '# purity_cval' in line:
                            if not isinstance(line.split('=')[1].strip(), float):
                                cur_cval = line.split('=')[1].strip()
                            else:
                                cur_cval = float(line.split('=')[1].strip())
                    else:
                        print("Error in ParseOut: hisens_vs_purity has invalid value. Must be 'hisens' or 'purity'.\nConfirm FacetsMeta has the value properly initialized.")
                        sys.exit()

                    #Get ploidy value.
                    if '# Ploidy' in line:
                        if not isinstance(line.split('=')[1].strip(), float):
                            cur_ploidy = line.split('=')[1].strip()
                        else:
                            cur_ploidy = float(line.split('=')[1].strip())

                    #Get dipLogR value.
                    if '# dipLogR' in line:
                        if not isinstance(line.split('=')[1].strip(), float):
                            cur_dipLogR = line.split('=')[1].strip()
                        else:
                            cur_dipLogR = float(line.split('=')[1].strip())

            return [cur_purity,cur_cval,cur_ploidy,cur_dipLogR]

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError parsing .out file. Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()



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
                cur_gene_level_data.append([gene_level_df['gene'][index], gene_level_df['gene_start'][index], gene_level_df
                                            ['gene_end'][index], gene_level_df['seg_start'][index], gene_level_df['seg_end'][index], 
                                            gene_level_df['seg_length'][index], gene_level_df['cf.em'][index], gene_level_df['tcn.em'][index], 
                                            gene_level_df['lcn.em'][index], gene_level_df['cn_state'][index], gene_level_df['filter'][index],
                                            gene_level_df['tsg'][index], gene_level_df['seg'][index], gene_level_df['median_cnlr_seg'][index],
                                            gene_level_df['segclust'][index], gene_level_df['mcn'][index], gene_level_df['genes_on_seg'][index],
                                            gene_level_df['gene_snps'][index], gene_level_df['gene_het_snps'][index], gene_level_df['spans_segs'][index],
                                            ])
            
            return cur_gene_level_data

        except Exception as e:
            print (bcolors.FAIL)
            print ("\tError parsing sample qc file. Terminating execution.")
            print ("\t" + gene_level_file)
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    ######################
    # parseSampleQC:  This function accepts a path to a .sample.qc.txt file and an expected cval
    #                 which can be identified with the parseOut function.  It returns a list
    #                 containing data for whole genome duplication, fraction genome altered, and fraction LOH.
    #                 Returns: [wgd,fga,frac_loh]
    ######################
    @staticmethod
    def parseSampleQC(qc_file, expected_cval):
        try:
            #Extract QC data.
            qc_df   = pd.read_csv(qc_file, sep="\t", low_memory=False)

            #Sometimes there are missing values, we need to handle those properly.
            #I.E. there is no fga column in refit_c50_pc100_diplogR_-0.25/P-0000543-T01-IM3_P-0000543-N01-IM3.qc.txt
            if 'fga' not in qc_df.columns or 'cval' not in qc_df.columns or 'frac_loh' not in qc_df.columns:
                return False
            else:
                #We want the row with a cval matching our target sample/fit.
                target_qc_row = qc_df.loc[qc_df['cval'] == int(expected_cval)]
                cur_wgd       = bool(target_qc_row['wgd'].values[0])
                cur_fga       = float(target_qc_row['fga'].values[0])
                cur_frac_loh  = float(target_qc_row['frac_loh'].values[0])

                return [cur_wgd,cur_fga,cur_frac_loh]

        except Exception as e:
            print (bcolors.FAIL)
            print ("\tError parsing sample qc file. Terminating execution.")
            print ("\t" + qc_file)
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    ######################
    # parseFacetsQC:  This function accepts the directory map list for a sample.
    #                 It will extract and return the "facets_qc" value for the target sample.
    ######################
    @staticmethod
    def parseFacetsQC(sample_dir_map):
        try:
            #Extract FACETS QC data.
            qc_df   = pd.read_csv(sample_dir_map[3], sep="\t", low_memory=False)

            #Sometimes there are missing values, we need to handle those properly.
            if 'fit_name' not in qc_df.columns or 'facets_qc' not in qc_df.columns:
                print("parseFacetsQC Fail")
                print(sample_dir_map)
                print(qc_df)
                return -1
            else:
                #Pick the row with the fit we've determined to be the best for the given sample.
                target_fit_name   = sample_dir_map[0].split("/")[-2]
                target_qc_row     = qc_df.loc[qc_df['fit_name'] == str(target_fit_name)]
                cur_pass_qc       = bool(target_qc_row['facets_qc'].values[0])

                return cur_pass_qc

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError parsing facets qc file. Terminating execution.")
            print("Items in sample_dir_map: " + str(len(sample_dir_map)))
            print (sample_dir_map)
            print (e)
            print (bcolors.ENDC)
            # sys.exit()


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

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError processing clinical data (parseCNCF). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()


    ######################
    # buildFacetsRun:  This function accepts a sample ID and parses all relevant files and assembles
    #                     all desired values.  It creates and returns an object of type FacetsRun. 
    ######################
    @staticmethod
    def buildFacetsRun(sample_id, facets_metadata):
        try:
            #print ("Building facets sample for " + sample_id)
            cur_dir_map            = facets_metadata.master_file_dict.get(sample_id)

            #If the file exists, but is empty, return a failure.
            # [out_file, cncf_file, qc_file, facets_qc_file, selected_fit_dir, gene_level_file]}
            if os.path.getsize(cur_dir_map[MetaDictMap.OUT_FILE]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " out_file is empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[MetaDictMap.CNCF_FILE]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " cncf_file is empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[MetaDictMap.QC_FILE]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " qc_file is empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[MetaDictMap.FACETS_QC_FILE]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " facets_qc_file is empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[MetaDictMap.GENE_FILE]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " gene-level file is empty!")
                    print (bcolors.ENDC)
                return False

            #Get data from our relevant input files.
            out_data               = FacetsRun.parseOut(cur_dir_map[MetaDictMap.OUT_FILE])
            sample_qc_data         = FacetsRun.parseSampleQC(cur_dir_map[MetaDictMap.QC_FILE], out_data[1])
            facets_qc_out          = FacetsRun.parseFacetsQC(cur_dir_map)

            #If we were missing information in any of our files, we don't want to proceed with this sample.
            if facets_qc_out == -1:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " facets_qc_file parse failed when building FacetsRun object!")
                    print ("\t\t\tDirMap for this sample: " + str(cur_dir_map))
                    print (bcolors.ENDC)
                FacetsMeta.failed_samples.append(sample_id)
                return False
            if sample_qc_data is False:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print(sample_id)
                    print ("\t\tWarning: " + sample_id + " sample_qc_data parse failed when building FacetsRun object!")
                    print ("\t\t\tDirMap for this sample: " + str(cur_dir_map))
                    print (bcolors.ENDC)
                FacetsMeta.failed_samples.append(sample_id)
                return False
                
            #Get data from our data mappings.
            short_id = ""
            if facets_metadata.selectSingleRun is False:
                short_id = sample_id.split("_")[0]
                cur_cancer_type        = facets_metadata.cancer_type_map.get(short_id)
                cur_cancer_type_detail = facets_metadata.cancer_type_detail_map.get(short_id)
                cur_clin_purity        = facets_metadata.clinical_purity_map.get(short_id)
                cur_onko_code          = facets_metadata.onkotree_code_map.get(short_id)
                cur_patient_id         = facets_metadata.patient_id_map.get(short_id)
                cur_tmb                = facets_metadata.cvr_tmb_score_map.get(short_id)
                cur_msi                = facets_metadata.msi_score_map.get(short_id)
                review_status          = cur_dir_map[MetaDictMap.FIT_STATUS].split("#")  #reviewstatus#fit_name#reviewed_by
                cur_review_status      = review_status[0]
                fit_name               = review_status[1]
                reviewed_by            = review_status[2]
            else:
                cur_cancer_type        = facets_metadata.cancer_type_map.get(sample_id)
                cur_cancer_type_detail = facets_metadata.cancer_type_detail_map.get(sample_id)
                cur_clin_purity        = facets_metadata.clinical_purity_map.get(sample_id)
                cur_onko_code          = facets_metadata.onkotree_code_map.get(sample_id)
                cur_patient_id         = facets_metadata.patient_id_map.get(sample_id)
                cur_tmb                = facets_metadata.cvr_tmb_score_map.get(sample_id)
                cur_msi                = facets_metadata.msi_score_map.get(sample_id)
                review_status          = cur_dir_map[MetaDictMap.FIT_STATUS].split("#")  #reviewstatus#fit_name#reviewed_by
                cur_review_status      = review_status[0]
                fit_name               = review_status[1]
                reviewed_by            = review_status[2]

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
                                            FacetsMeta.hisens_vs_purity,
                                            cur_review_status,
                                            fit_name,
                                            reviewed_by
                                            )

            #Now build in the segments from the CNCF file.
            segments = FacetsRun.parseCNCF(cur_dir_map[MetaDictMap.CNCF_FILE],cur_dir_map[MetaDictMap.ADJUSTED_SEG_FILE])
            for i in range(len(segments)):
                if str(segments[i].cf) == "nan" or str(segments[i].cf_base) == "nan":
                    continue
                cur_facets_sample.addSegment(segments[i])

            #Now build in the gene-level data.
            geneLevel = FacetsRun.parseGeneLevel(cur_dir_map[MetaDictMap.GENE_FILE])
            for curGene in geneLevel:
                geneToAdd = FacetsGene(curGene[0],
                                        curGene[1],
                                        curGene[2],
                                        curGene[3],
                                        curGene[4],
                                        curGene[5],
                                        curGene[6],
                                        curGene[7],
                                        curGene[8],
                                        curGene[9],
                                        ###
                                        curGene[10],
                                        curGene[11],
                                        curGene[12],
                                        curGene[13],
                                        curGene[14],
                                        curGene[15],
                                        curGene[16],
                                        curGene[17],
                                        curGene[18],
                                        curGene[19],
                                        )

                cur_facets_sample.addGene(geneToAdd)

            return cur_facets_sample

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.buildFacetsRun(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

        


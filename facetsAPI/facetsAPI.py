from heapq import merge
from random import sample
import sys
import os.path
import glob
import pandas as pd
import pickle


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
# FacetsMeta:    This class represents necessary metadata structures that hold information such as directory paths,
#                cancer types, and other clinical data that need to be referenced in FacetsRuns.
######################
class FacetsMeta:
    use_unreviewed_defaults = False
    hisens_vs_purity        = "purity"
    selectSingleRun         = False
    failed_samples          = []
    run_verbose             = False

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

    def __init__(self, clinical_sample_file, facets_repo_path, hisens_vs_purity, persist_data="no"):
        #Relevant path and file data.
        self.clinical_sample_file    = clinical_sample_file
        self.facets_repo_path        = facets_repo_path
        self.hisens_vs_purity        = hisens_vs_purity
        self.persist_data            = persist_data

        #Data structures and storage.
        self.master_file_dict       = {} # A map of relevant files for each sample. {id: [out_file, cncf_file, qc_file, facets_qc_file, selected_fit_dir]}

        #These come from data_clinical_sample.
        self.cancer_type_map        = {} # A map of sample ids to cancer types.
        self.cancer_type_detail_map = {} # A map of sample ids to cancer type detailed.
        self.patient_id_map         = {} # A map of sample ids to patient ids.
        self.clinical_purity_map    = {} # A map of purity values from the clinical sample file.
        self.onkotree_code_map      = {} # A map of onkotree codes.
        self.cvr_tmb_score_map      = {} # A map of cvr tmb scores from the clinical sample file.
        self.msi_score_map          = {} # A map of msi scores.
        self.long_id_map            = {} # A map of sample ids to their corresponding long_ids.  id -> [long_id1, long_id2...]

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

    #This function will set this FacetsMeta object to select a single run per sample.
    def setSingleRunPerSample(self, doSet):
        try:
            if not isinstance(doSet, bool):
                print (bcolors.FAIL)
                print ("\t\tError in FacetsMeta.setSingleRunPerSample(). Expected True or False value.")
                print (e)
                print (bcolors.ENDC)
                sys.exit()
            FacetsMeta.selectSingleRun = doSet
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsMeta.setSingleRunPerSample(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

    #This function will set this FacetsMeta object to allow default fits when selecting single runs.
    def allowDefaultFitsIfNoBest(self, doSet):
        try:
            if not isinstance(doSet, bool):
                print (bcolors.FAIL)
                print ("\t\tError in FacetsMeta.allowDefaultFitsIfNoBest(). Expected True or False value.")
                print (e)
                print (bcolors.ENDC)
                sys.exit()
            if not FacetsMeta.selectSingleRun and doSet:
                print (bcolors.FAIL)
                print("\t\tError in FacetsMeta.allowDefaultFitsIfNoBest().  Cannot set default fit selection to True when FacetsMeta.setSingleRunPerSample() is False.\nRun FacetsMeta.setSingleRunPerSample(True) first.")
                print("\t\t\tFacetsMeta.selectSingleRun: " + str(FacetsMeta.selectSingleRun))
                print("\t\t\tFailed to set FacetsMeta.use_unreviewed_defaults to: " + str(FacetsMeta.doSet))
                print (bcolors.ENDC)
                sys.exit()

            FacetsMeta.use_unreviewed_defaults = doSet
        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsMeta.allowDefaultFitsIfNoBest(). Terminating execution.")
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
            for id in target_ids:
                cur_short_id = id[0:7]
                cur_id_dirs = glob.glob(self.facets_repo_path + cur_short_id + "/" + id + "*")

                #Keep track of how many samples had the same tumor/id, but multiple normals.
                if len(cur_id_dirs) > 1:
                    num_multi_normals = num_multi_normals + 1

                if not cur_id_dirs:
                    continue
                else:
                    id_with_normal = cur_id_dirs[-1].split('/')[-1]

                    cur_manifest_file  = cur_id_dirs[-1] + "/facets_review.manifest"
                    cur_facets_qc_file = cur_id_dirs[-1] + "/facets_qc.txt" 

                    #If there is no manifest file, skip the sample.
                    if not os.path.exists(cur_manifest_file):
                        num_missing_manifest = num_missing_manifest + 1
                        continue

                    cur_manifest = pd.read_csv(cur_manifest_file, skiprows=1, sep="\t", low_memory=False)
                    manifest_data = cur_manifest[['path', 'review_status', 'fit_name']]

                    #If we are looking for individual fits we will want to look for the best and consider default fits.
                    if FacetsMeta.selectSingleRun:
                        #Grab important data for best/acceptable/default fits.
                        target_best_fit = manifest_data.loc[manifest_data['review_status'].isin(["reviewed_best_fit"])]
                        target_accept_fit = manifest_data.loc[manifest_data['review_status'].isin(["reviewed_acceptable_fit"])]
                        if(FacetsMeta.use_unreviewed_defaults):
                            target_default_fit = manifest_data.loc[manifest_data['fit_name'].isin(["default"])]

                        selected_fit_dir = ""
                        #If there is a best fit available, we want to use that.
                        if(not target_best_fit.empty):
                            best_fit_list = target_best_fit.values.tolist()[0]
                            if best_fit_list[0][-1] != '/':
                                best_fit_list[0] = best_fit_list[0] + "/"
                            selected_fit_dir  = best_fit_list[0] + best_fit_list[2] + "/"
                            num_best_fits = num_best_fits + 1     
                        #If there is no best fit, check for acceptable fits.
                        elif(not target_accept_fit.empty):
                            accept_fit_list = target_accept_fit.values.tolist()[0]
                            if accept_fit_list[0][-1] != '/':
                                accept_fit_list[0] = accept_fit_list[0] + "/"
                            selected_fit_dir  = accept_fit_list[0] + accept_fit_list[2] + "/"
                            num_acceptable_fits = num_acceptable_fits + 1
                        #If we have nothing still, and are accepting default fits, use that.
                        elif(FacetsMeta.use_unreviewed_defaults):
                            default_fit_list = target_default_fit.values.tolist()[0]
                            if default_fit_list[0][-1] != '/':
                                default_fit_list[0] = default_fit_list[0] + "/"
                            selected_fit_dir  = default_fit_list[0] + default_fit_list[2] + "/"
                            num_default_fits = num_default_fits + 1
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
                        if self.hisens_vs_purity == "purity":
                            cur_out    = glob.glob(selected_fit_dir + "/*_purity.out")
                            cur_cncf   = glob.glob(selected_fit_dir + "/*_purity.cncf.txt")
                        elif self.hisens_vs_purity == "hisens":
                            cur_out    = glob.glob(selected_fit_dir + "/*_hisens.out")
                            cur_cncf   = glob.glob(selected_fit_dir + "/*_hisens.cncf.txt")
                        else:
                            print("Error: hisens_vs_purity value should be 'hisens' or 'purity'.")
                            sys.exit()

                        cur_qc         = glob.glob(selected_fit_dir + "*qc.txt")
                        cur_gene_level = glob.glob(selected_fit_dir + "*gene_level.txt")

                        #If critical files are missing, skip the sample.
                        if not cur_out or not cur_cncf or not cur_qc or not cur_gene_level:
                            num_missing = num_missing + 1
                            continue
                        else:
                            out_file        = cur_out[0]
                            cncf_file       = cur_cncf[0]
                            qc_file         = cur_qc[0]
                            gene_level_file = cur_gene_level[0]

                        self.long_id_map[id]   = [id_with_normal]
                        self.master_file_dict[id] = [out_file, cncf_file, qc_file, cur_facets_qc_file, selected_fit_dir, gene_level_file]

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

                            #Handle cases where the path doesn't properly end with a / character.
                            cur_sample_folder = row['path']
                            if cur_sample_folder[-1] != '/':
                                cur_sample_folder = cur_sample_folder + "/"

                            cur_fit_folder = cur_sample_folder + row['fit_name'] + "/"

                            if self.hisens_vs_purity == "purity":
                                cur_out    = glob.glob(cur_fit_folder + "*_purity.out")
                                cur_cncf   = glob.glob(cur_fit_folder + "*_purity.cncf.txt")
                            elif self.hisens_vs_purity == "hisens":
                                cur_out    = glob.glob(cur_fit_folder + "*_hisens.out")
                                cur_cncf   = glob.glob(cur_fit_folder + "*_hisens.cncf.txt")
                            else:
                                print("Error: hisens_vs_purity value should be 'hisens' or 'purity'.")
                                sys.exit()

                            cur_qc         = glob.glob(cur_fit_folder + "*qc.txt")
                            cur_gene_level = glob.glob(cur_fit_folder + "*gene_level.txt")

                            #If critical files are missing, skip the sample.
                            if not cur_out or not cur_cncf or not cur_qc or not cur_gene_level:
                                num_missing = num_missing + 1
                                continue
                            else:
                                out_file        = cur_out[0]
                                cncf_file       = cur_cncf[0]
                                qc_file         = cur_qc[0]
                                gene_level_file = cur_gene_level[0]

                            cur_run_list.append(long_id)
                            self.master_file_dict[long_id] = [out_file, cncf_file, qc_file, cur_facets_qc_file, cur_fit_folder, gene_level_file]

                        self.long_id_map[id] = cur_run_list

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
            print ("\t\tError processing clinical data. Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()



######################
# FacetsDataset:    This class represents a set of FacetsRuns.  It includes functions for selecting, filtering,
#                   writing output, and modifying/creating file structures.
######################
class FacetsDataset:
    def __init__(self):
        self.sampleList = {}
        self.runList = []
        self.numSamples = 0
        self.ref_facetsMeta = None

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

    #This function will write a facets data set to a tab delimited file.
    def writeToDelimFile(self, include_):
        pass

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
                    print("hi" + curSample.cancerType)
                    curSample.printSample()
                    curLineString = ""
                    curLineString += curSample.id + "\t"
                    curLineString += curSample.cancerType + "\t"
                    curLineString += curSample.cancerTypeDetail + "\t"
                    curLineString += str(curSample.purity) + "\t"
                    curLineString += curSample.onkoCode + "\t"
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
    def buildFacetsDataset(self, facets_metadata):
        try:
            if not isinstance(facets_metadata, FacetsMeta):
                print (bcolors.FAIL)
                print("Error in FacetsDataset.buildFacetsDataset(): buildFacetsDataset must be an initialized object of type FacetsMeta.")
                print (bcolors.ENDC)
                sys.exit()

            self.ref_facetsMeta = facets_metadata

            total_samples_prepped = 0
            num_build_failed      = 0
            num_file_failed       = 0
            print("\tApplying filtering and building FacetsDataset beginning with " + str(len(facets_metadata.master_file_dict)) + " samples.")

            #Go through our master file dictionary and build facets sample objects.
            for key in facets_metadata.master_file_dict.keys():
                #print(key)
                if total_samples_prepped % 1000 == 0 and total_samples_prepped != 0:
                    print("\t\tSamples Processed: " + str(total_samples_prepped))

                cur_dir_map = facets_metadata.master_file_dict.get(key)

                #CancerType Filter.
                if self.filterCancerType:
                    curCancerType = facets_metadata.cancer_type_map.get(key)
                    if curCancerType not in self.selectedCancerTypes:
                        continue

                #CancerTypeDetail Filter.
                if self.filterCancerTypeDetail:
                    curCancerDetailType = facets_metadata.cancer_type_detail_map.get(key)
                    if curCancerDetailType not in self.curCancerDetailType:
                        continue

                #Clinical Purity Filter.
                if self.filterClinicalPurity:
                    curClinPurity = facets_metadata.clinical_purity_map.get(key)
                    if curClinPurity == "NA":
                        continue
                    curClinPurity = float(curClinPurity)
                    if curClinPurity < self.minClinPurity or curClinPurity > self.maxClinPurity:
                        continue

                #OnkoCode Filter.
                if self.filterOnkoCode:
                    curOnkoCode = facets_metadata.onkotree_code_map.get(key)
                    if curOnkoCode not in self.selectedOnkoCodes:
                        continue

                #TMB Filter.
                if self.filterTMB:
                    cur_tmb = facets_metadata.cvr_tmb_score_map.get(key)
                    if cur_tmb == "NA":
                        continue
                    cur_tmb = float(cur_tmb)
                    if cur_tmb < self.minTMB or cur_tmb > self.maxTMB:
                        continue

                #MSI Filter.
                if self.filterMSI:
                    cur_msi = facets_metadata.msi_score_map.get(key)
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
                cur_run = FacetsRun.buildFacetsRun(key,facets_metadata)
                if cur_run is False:
                    num_build_failed += 1
                    continue
                else:
                    cur_sample_id = ""
                    if FacetsMeta.selectSingleRun:
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

                                cur_sample.printSample()
                        else:
                            print (bcolors.FAIL)
                            print ("\t\tError in FacetsDataset.buildFacetsDataset(). No long id found for " + str(key))
                            print (e)
                            print (bcolors.ENDC)
                            sys.exit()
                    else:
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

    #This function sets a facets purity filter.  The selectedCancerDetailTypes parameter should be a list of strings.
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
    def __init__(self, gene, gene_start, gene_end, seg_start, seg_end, seg_length, cf, tcn, lcn, cn_state, filter):
        self.gene       = gene
        self.gene_start = gene_start
        self.gene_end   = gene_end
        self.seg_start  = seg_start
        self.seg_end    = seg_end
        self.seg_length = seg_length
        self.cf         = cf
        self.tcn        = tcn
        self.lcn        = lcn
        self.cn_state   = cn_state
        self.filter     = filter

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

    def __init__(self, chrom, start, end, cnlr_median, cf, cf_base, tcn, lcn):
        # These values are calculated in FacetsRun.defineArms().
        self.arm         = "undefined" # This is the chr/arm.  I.E. 1p, 5q, etc.
        self.percentArm  = -1          # This is the percentage of the arm that this segment covers. 

        self.chrom       = int(chrom)
        self.start       = int(start)
        self.end         = int(end)
        self.cnlr_median = float(cnlr_median)
        self.length      = int(max(end - start, start - end))
        self.cf          = float(cf)      #This API uses cf.em values for everything, so its just called cf here.
        
        try:
            self.cf_base = float(cf_base)
        except ValueError:
            try:
                self.cf_base = int(cf_base)
            except ValueError:
                self.cf_base = -1

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
    

    def printSample(self):
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)
        print("FacetsSample: " + str(self.id))
        print("\tThis sample contains " + str(len(self.runs)) + " FacetsRuns.")
        for i in range(len(self.runs)):
            print("ID: " + str(self.runs[i].id))
            print("Fit Dir: " + str(self.runs[i].fitDir))
        print(bcolors.BOLD + "\t~-===--===--===--===--===--===--===--===-~" + bcolors.ENDC)


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
                 wgd, fga, frac_loh, facets_qc, tmb, msi, purity_hisens):
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
                        cur_seg.percentArm = float(cur_seg.length) / float(cur_chrom_ranges[3])
                    #If not, split the segment.
                    else:
                        segsToRemove.append(cur_seg)
                        
                        #Split into two segments around the centromere.
                        seg1_newStart = cur_seg.start
                        seg1_newEnd   = cur_chrom_ranges[1] - 1
                        new_Seg1      = FacetsSegment(cur_seg.chrom,seg1_newStart,seg1_newEnd,cur_seg.cnlr_median,cur_seg.cf,cur_seg.cf_base,cur_seg.tcn,cur_seg.lcn)
                        new_Seg1.arm  = str(cur_seg.chrom) + "p"
                        new_Seg1.percentArm = float(new_Seg1.length) / float(cur_chrom_ranges[1])

                        seg2_newStart = cur_chrom_ranges[2] + 1
                        seg2_newEnd   = cur_seg.end
                        new_Seg2      = FacetsSegment(cur_seg.chrom,seg2_newStart,seg2_newEnd,cur_seg.cnlr_median,cur_seg.cf,cur_seg.cf_base,cur_seg.tcn,cur_seg.lcn)
                        new_Seg2.arm  = str(cur_seg.chrom) + "q"
                        new_Seg2.percentArm = float(new_Seg2.length) / float(cur_chrom_ranges[3])

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
                                            gene_level_df['lcn.em'][index], gene_level_df['cn_state'][index], gene_level_df['filter'][index]])
            
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
            sys.exit()


    ######################
    # parseCNCF:      This function accepts a CNCF file path and parses it for data.
    #                 It will build and return a list of FacetsSegment type objects representing the CNCF.
    ######################
    @staticmethod
    def parseCNCF(cncf_file):
        try:
            cncf_df = pd.read_csv(cncf_file, sep="\t", low_memory=False)

            #Extract the relevant columns.
            chrom       = cncf_df['chrom'].tolist()
            loc_start   = cncf_df['loc.start'].tolist()
            loc_end     = cncf_df['loc.end'].tolist()
            clnr_median = cncf_df['cnlr.median'].tolist()
            cf_em       = cncf_df['cf.em'].tolist()
            cf_base     = cncf_df['cf'].tolist()
            tcn_em      = cncf_df['tcn.em'].tolist()
            lcn_em      = cncf_df['lcn.em'].tolist()

            #Build a FacetsSegment object for each row.
            seg_list = []
            for i in range(len(chrom)):
                cur_seg = FacetsSegment(chrom[i],
                                        loc_start[i],
                                        loc_end[i],
                                        clnr_median[i],
                                        cf_em[i],
                                        cf_base[i],
                                        tcn_em[i],
                                        lcn_em[i])
                seg_list.append(cur_seg)

            return seg_list

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError processing clinical data. Terminating execution.")
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
            # [out_file, cncf_file, qc_file, facets_qc_file, selected_fit_dir]}
            if os.path.getsize(cur_dir_map[0]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " out_file is 0 empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[1]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " cncf_file is 0 empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[2]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " qc_file is 0 empty!")
                    print (bcolors.ENDC)
                return False
            if os.path.getsize(cur_dir_map[3]) == 0:
                if facets_metadata.run_verbose:
                    print (bcolors.WARNING)
                    print ("\t\tWarning: " + sample_id + " facets_qc_file is 0 empty!")
                    print (bcolors.ENDC)
                return False

            #Get data from our relevant input files.
            out_data               = FacetsRun.parseOut(cur_dir_map[0])
            sample_qc_data         = FacetsRun.parseSampleQC(cur_dir_map[2], out_data[1])
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
            else:
                cur_cancer_type        = facets_metadata.cancer_type_map.get(sample_id)
                cur_cancer_type_detail = facets_metadata.cancer_type_detail_map.get(sample_id)
                cur_clin_purity        = facets_metadata.clinical_purity_map.get(sample_id)
                cur_onko_code          = facets_metadata.onkotree_code_map.get(sample_id)
                cur_patient_id         = facets_metadata.patient_id_map.get(sample_id)
                cur_tmb                = facets_metadata.cvr_tmb_score_map.get(sample_id)
                cur_msi                = facets_metadata.msi_score_map.get(sample_id)

            #Build the FacetsRun object.
            cur_facets_sample = FacetsRun(sample_id,
                                            cur_patient_id, 
                                            cur_dir_map[4], 
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
                                            FacetsMeta.hisens_vs_purity)
            
            #Now build in the segments from the CNCF file.
            segments = FacetsRun.parseCNCF(cur_dir_map[1])
            for i in range(len(segments)):
                if str(segments[i].cf) == "nan" or str(segments[i].cf_base) == "nan":
                    continue
                cur_facets_sample.addSegment(segments[i])

            return cur_facets_sample

        except Exception as e:
            print (bcolors.FAIL)
            print ("\t\tError in FacetsRun.buildFacetsRun(). Terminating execution.")
            print (e)
            print (bcolors.ENDC)
            sys.exit()

        


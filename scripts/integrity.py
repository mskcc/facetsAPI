import os

import pandas as pd
import sys
import os
from datetime import datetime
#change this to wherever the facetsAPI is stored
# sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')
sys.path.insert(1, '/juno/work/ccs/orgeraj/facetsAPI')

from facetsAPI import *

def integrity_check(inFacetsMeta = False, inFacetsDataset = False):
    """
    This script will check the integrity of the entered Facets dataset is it has the proper file structure.
        Any empty files (size 0 bytes). 
        Any missing files (just don't exist)
        All files have the correct number of columns. (For example, some files qc file were missing WGD columns in rare cases).-
        All files have existing bam files that in an expected size range. (Some bams run on aug 9 are corrupted and failed in a recent run, although they worked in previous runs.)-
        Manifest files do not have multiple fits of the same type. (I.e. should not have 2+ "best fit" options)
        All fits listed in the manifest file should have corresponding fit folders, and associated files.
        All files listed in the cohort_level manifest file have an associated folder.
        Sample still has consent. (12245-a)
        Number of fits
        Samples that exist in the repository, but are not listed in data_clinical_sample.

    Metrics
        Total number of samples
        Total number of fits.
        Number of best/acceptable/default fits.
        Number of samples with multiple normals but the same tumor.

    Inputs: 
    inFacetsMeta (FacetsMeta object):  Inputted facets meta object that has been constructed previously.  If no 
                                        inFacetsDataset is input then the integrity check will check the facets directory
                                        indicated in FacetsMeta
    inFacetsDataset (facetsDataset object) : Inputted facetsDataset object that has been constructed previously. The value
                                                inFacetsDataset.ref_facetsMeta will replace te input FacetsMeta object if there 
                                                was one. 
    Outputs:
    FACETSintegrity_{datetime}.tsv : A TSV containing the data collected from integrity.py
    Cols are as such : 
    ID : Name of the sample folder
    emptyfiles:  Files that have a size of 0 bytes but do exist
    missing_filesin_samplefolder: Files that do not exist in the sample folder and should be there
    MissingfilesinFitfolder:  Files that do not exist in the fit folder but should be there
    missingcols:  Files that are missing columns from the fit folder
    tumorbam_size_pass: The tumor bam file is within two standard deviations of the mean file size
    normalbam_size_pass: The tumor bam file is within two standard deviations of the mean file size
    manifest_file_labelerror: Manifest files that have multiple fits of the same type
    manifest_folder_missing:  A fit in the manifest does not have the fit folder
    Consent:  12245-a consent.  1 if consented
    FitCount: Number of fits
    In_Clinical_file:  If constructed from facetsMeta this will be 1 if the sample is in the clincal file. 

    Also in terminal are printedseveral metrics such as: 
    Total number of samples, 
    number of fits, 
    Number of Best, acceptable, default, and no fit
    Number of samples with multiple normals but same tumor                                                  
    """ 


    EvalDirectory = False
    if inFacetsDataset :
        inFacetsMeta = inFacetsDataset.ref_facetsMeta
    elif inFacetsMeta:
        EvalDirectory = True
    else:
        print (bcolors.FAIL)
        print ("\t\tError in integrity_check. Terminating execution. Either FacetsMeta must be provided or Facets Dataset")
        print (e)
        print (bcolors.ENDC)
        sys.exit()


    patientIDorganized = {}
    sample_num = 0
    fit_num=0   
    fitstats=[0,0,0,0,0] #best fit,acceptable,default fit, not reviewed, no fit
    manyNormal_sameTumor = 0
    normaldict = {}  # Checks for manyNormal_sameTumor
    startdir = inFacetsMeta.facets_repo_path



    #create hashtable of bams connected to items 
    keyfile = "/juno/work/ccs/shared/resources/impact/dmp_mirror_key/dmp_key_paired_latest.txt"
    keydict={}
    header = True
    print("Creating sizedict")
    with open(keyfile,"r") as keyf: 
        
        for line in keyf:
            if header:
                header =False
                pass
            else:
                    
                line = line.split("\t")
                # print(line[4])
                
                # print(line[10])

                if line[0] not in keydict and line[4] != "NA":
                    
                    path = "/juno/res/dmpcollab/dmpshare/share/irb12_245/" + line[4][0] + "/" + line[4][1] + "/" + line[4] + ".bam"
                    keydict[line[0]] = [line[4],path]  #[filename,path]

                if line[8] not in keydict and line[10] != "NA": 
                    path = "/juno/res/dmpcollab/dmpshare/share/irb12_245/" + line[10][0] + "/" + line[10][1] + "/" + line[10] + ".bam"
                    keydict[line[9]] = [line[10],path]
    print("Sizedict done")

    #create consent hashtable
    print("creating consent table")
    consent_dict = {}
    consenttable = "/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/dmp-2022/mskimpact/data_clinical_patient.txt"
    with open(consenttable,"r") as conf: 
        header = 0
        for line in conf:
            if header<5:
                header +=1
                
            else:
                line = line.split("\t")
                consent_dict[line[0]] = line[7]
                print(line[7])
    print("consent done")
    #######################################


    # Create data_clincal_sample hashtable to check if patientID is in data_clinical_sample
    if inFacetsMeta.hasClinicalMeta:
        print("data_clincial table creation")
        # clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
        header = True
        clinical_sample_dict = {}
        with open(inFacetsMeta.clinical_sample_file,"r") as clinf:
            if header:
                    header =False
                    pass
            else:
                line = line.split("\t")
                clinical_sample_dict.add(line[0]) 
                
        print("data_clincial table done")
        ############################################


    if EvalDirectory: 
        startdir = os.listdir(inFacetsMeta.facets_repo_path)
        print(startdir[0])
        searchdirectories = [inFacetsMeta.facets_repo_path+x for x in startdir]
        print(searchdirectories[0])
        fullsampledirectories = []
        for item in searchdirectories:
            if os.path.isdir(item):
                tmpli = os.listdir(item)
                # print(tmpli )
                for k in tmpli:
                    fullsampledirectories.append(item+"/"+k+"/")

                # [fullsampledirectories.append(item+k) for k in tmpli]


        
    else:
        startdir=""
        fullsampledirectories = [inFacetsMeta.facets_repo_path+"/"+x[0:7]+"/"+x+"/" for x in inFacetsDataset.sampleList]

    #####Start checks##############
    print(EvalDirectory)

    for samplefolder in fullsampledirectories:
        if "backup" in samplefolder or "test" in samplefolder:
            continue 
        # try:
        if os.path.isdir(samplefolder):
            if len(startdir) ==0:
                pass
            else: 
                pass    
            must_have_outerfiles = {"facets_qc.txt", "facets_review.manifest", ".maf", "nonsignedout.maf"}

            # Move into a individual sample folder 
            samplepath = samplefolder
            sample_num+=1  #add to count
            normalID = samplefolder.split("/")[-2].split("_")[1]
            tumorID = samplefolder.split("/")[-2].split("_")[0]    
            patientID = tumorID[:-8]
            print(tumorID)
            
            if samplefolder not in patientIDorganized:
                #Establish the dictionary
                patientIDorganized[samplefolder] = [[], #empty files
                                                    [], #Missing files from samplefolder
                                                    [], #Missing files from fit folder
                                                    [], #Files missing cols
                                                    [0,0], #Tumorbam size, Normalbam size. 1 is passing, accpetable size
                                                    [], #1 if multiple best fits found
                                                    [], #Missing fit folders
                                                    0, #Has consent: 1 if consented
                                                    0, #number of fits
                                                    [], #fits with two entries in Manifest
                                                    0] #in clinical folder: 1 if so

            if normalID not in normaldict: 
                normaldict[normalID] = 1
            else:
                manyNormal_sameTumor +=1
                normaldict[normalID] +=1


            # patientIDorganized[samplefolder][0] = samplefolder
            if os.path.isfile(samplefolder+"/facets_review.manifest"):
                must_have_outerfiles.remove("facets_review.manifest")

                manifest =  pd.read_csv(samplefolder+"/facets_review.manifest",sep="\t", header = 1)

                fit_name_li =set()
                # fit_name_set = []
                best_fit_found = False
                for index,row in manifest.iterrows():
                    if row["fit_name"] != "Not selected":
                        if row["fit_name"] not in fit_name_li:
                            fit_name_li.add(row["fit_name"])
                        else:
                            patientIDorganized[samplefolder][9].append(row["path"]) #add something
                        
                        #checks if manifest folder contains multiple best fits
                        if best_fit_found and "best_fit" in row['review_status']:
                            patientIDorganized[samplefolder][5] = 1 #1 if multiple best fits
                        else:
                            best_fit_found ==True

                    #fitstats 
                    if "best" in row['review_status']:
                        fitstats[0]+=1
                    elif "acceptable" in row['review_status']:
                        fitstats[1]+=1
                    elif "default" in row['review_status']:
                        fitstats[2]+=1
                    elif "not_reviewed" in row['review_status']:
                        fitstats[3]+=1
                    elif "no_fit" in row['review_status']:
                        fitstats[4]+=1
                    # print(row['review_status'])

                    #Adjust fitnumber
                    patientIDorganized[samplefolder][8] +=1
            else:
                print("No Manifest file")

            #Check size of bam 1 means fail
            try:
                normsize = os.path.getsize(keydict[normalID][1])
                if normsize > 956598954: 
                    patientIDorganized[samplefolder][4][0] = 1
            except:
                print("Normal ID not found in dmp key")
                patientIDorganized[samplefolder][4][0] = "Could not find"


            try:
                tumorsize = os.path.getsize(keydict[tumorID][1])
                if tumorsize > 956598954: 
                    patientIDorganized[samplefolder][4][1] = 1
            except:
                print("Tumor ID not found in dmp key")
                patientIDorganized[samplefolder][4][1] = "Could not find"
            ##################
                
            #Check for consent in consent dict: 
            try:    
                if "YES" in consent_dict[patientID]:
                    patientIDorganized[samplefolder][7] = 1

            except: 
                patientIDorganized[samplefolder][7] = "Unknown"



            # Check if in data_clinical_sample

            if EvalDirectory:
                    
                if tumorID in clinical_sample_dict:
                    patientIDorganized[samplefolder][9] = 1

            else:
                patientIDorganized[samplefolder][9] = None

            try:
                for runfi in os.listdir(samplepath):
                    #Iterates through items in the sample folder.  We will use this to check for missing files
                    # as well as missing manifest folders, and eventually checking each individual file's cols
                    # print(("samplepath+runfi",samplepath,runfi))
                    #Establish Must have file list.  Items are removed from this list as they are found
                    must_have_files = {"hisens_diplogR.adjusted.seg", "hisens_diplogR.unadjusted.seg", "hisens.CNCF.png", "hisens.cncf.txt", "hisens.out","hisens.Rdata", 
                        "purity_diplogR.adjusted.seg", "purity_diplogR.unadjusted.seg", "purity.CNCF.png", "purity.cncf.txt", "purity.out","purity.Rdata", 
                        "arm_level.txt", ".ccf.maf","gene_level.txt",".nonsignedout.ccf.maf", ".qc.txt", "purity.seg", "hisens.seg"}

                    if os.path.isdir(samplepath+"/"+runfi):
                        #This moves into fit folders and checks each file for size and if any are missing
                        if runfi in fit_name_li:
                            fit_name_li.remove(runfi)

                        
                        for fi in os.listdir(samplepath+"/"+runfi):
                            if fi[0] == ".":
                                continue
                            #Moving into individuual fit folders

                            #Check filesize
                            if os.path.getsize(samplepath+"/"+runfi+"/"+fi) == 0 :
                                #add to empty files
                                patientIDorganized[samplefolder][0].append(samplepath+"/"+runfi+"/"+fi)

                            fi_li = fi.split('_')
                            print(fi_li)

                            if len(fi_li)>3:
                                #this takes care of hisens and purity_diplogR.adj style files
                                print(">3")
                                if fi.endswith("png"):
                                    try:
                                        must_have_files.remove(fi_li[2])
                                    except: 
                                        must_have_files.remove(fi_li[3])
                                else:
                                        
                                    must_have_files.remove(fi_li[2]+"_"+fi_li[3])
                                    with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                        for line in segfile:
                                            header=set(line.replace("\n","").split("\t"))
                                            break
                                    headtst ={"ID", "chrom", "loc.start", "loc.end", "num.mark","seg.mean"}
                                    if header==headtst:
                                        pass
                                    else:
                                        if len(header)> len(headtst):
                                            print(samplepath+"/"+runfi+"/"+fi)
                                        patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi)  #Files missing cols 
                            elif len(fi_li)==3:
                                print("==3")
                                if "level" in fi_li[2]:
                                    #gene_level
                                    # print(fi_li[1].split(".")[1]+"_"+fi_li[2])
                                    must_have_files.remove(fi_li[1].split(".")[1]+"_"+fi_li[2])

                                    if "gene" in fi_li[1].split(".")[1]:
                                        #Gene_level check
                                        with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                            for line in segfile:
                                                header=set(line.replace("\n","").split("\t"))
                                                break
                                        headtst = {"sample", "gene", "chrom", "gene_start", "gene_end", "tsg", "seg", "median_cnlr_seg", "segclust", "seg_start", "seg_end", 'cf.em',\
                                                "tcn.em", "lcn.em", "cf", "tcn", "lcn", "seg_length", "mcn", "genes_on_seg", "gene_snps", "gene_het_snps", "spans_segs", "cn_state", "filter"}
                                        if header==headtst:
                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 

                                    elif  "arm" in fi_li[1].split(".")[1]:
                                        #Arm_level check
                                        with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                            for line in segfile:
                                                header=set(line.replace("\n","").split("\t"))
                                                break
                                        headtst = {"sample", "arm", "tcn", "lcn", "cn_length", "arm_length", "frac_of_arm", "cn_state"}
                                        if header==headtst:
                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                                print(header)
                                                break
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 



                                else:
                                    # print("?")
                                    # print(fi_li[2])
                                    try:
                                        must_have_files.remove(fi_li[2])
                                    except:
                                        print("Extra file")
                                        pass

                                    # print(fi_li)

                                    if fi_li[2].endswith("out") or fi_li[2].endswith("png") or fi_li[2].endswith("Rdata"):
                                        pass
                                    else:
                                        with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                            for line in segfile:
                                                header=set(line.replace("\n","").split("\t"))
                                                # print(line.split("\t"))
                                                break
                                        headtst = {"ID", "chrom", "loc.start", "loc.end", "seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf", "tcn", "lcn", "cf.em", "tcn.em", "lcn.em"}
                                        if header==headtst:
                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 

                            if len(fi_li) == 2:
                                if fi.endswith(".qc.txt"):
                                    must_have_files.remove(".qc.txt")
                                    #.qc check
                                    with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                        for line in segfile:
                                            header=set(line.replace("\n","").split("\t"))
                                            break
                                        headtst ={"sample", "cval", "dipLogR_flag", "n_alternative_dipLogR", "wgd", "fga", "n_dip_bal_segs", "frac_dip_bal_segs", "n_dip_imbal_segs", "frac_dip_imbal_segs", "n_amps", \
                                                    "n_homdels", "frac_homdels", "n_homdels_clonal", "frac_homdels_clonal", "n_cn_states", "n_segs", "n_cnlr_clusters", "n_lcn_na", "n_loh", "frac_loh", "n_segs_subclonal", \
                                                    "frac_segs_subclonal", "n_segs_below_dipLogR", "frac_below_dipLogR", "n_segs_balanced_odd_tcn", "frac_balanced_odd_tcn", "n_segs_imbalanced_diploid_cn", \
                                                    "frac_imbalanced_diploid_cn", "n_segs_lcn_greater_mcn", "frac_lcn_greater_mcn", "n_snps" , "n_het_snps", "frac_het_snps", "n_snps_with_300x_in_tumor", \
                                                    "n_het_snps_with_300x_in_tumor", "n_het_snps_hom_in_tumor_1pct", "n_het_snps_hom_in_tumor_5pct", "frac_het_snps_hom_in_tumor_1pct", "frac_het_snps_hom_in_tumor_5pct", \
                                                    "mean_cnlr_residual", "sd_cnlr_residual", "n_segs_discordant_tcn", "frac_discordant_tcn", "n_segs_discordant_lcn", "frac_discordant_lcn", "n_segs_discordant_both", \
                                                    "frac_discordant_both", "n_segs_icn_cnlor_discordant", "frac_icn_cnlor_discordant", "mafr_median_all", "mafr_median_clonal", "mafr_n_gt_1"}
                                        if header== headtst:

                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                                print(header.difference(headtst))
                                                print("qc.txt")
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 
                                    
                                elif fi.endswith(".nonsignedout.ccf.maf"):
                                    must_have_files.remove(".nonsignedout.ccf.maf")
                                    # .nonsignedout check
                                    with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                        for line in segfile:
                                            header=set(line.replace("\n","").split("\t"))
                                            break
                                        headtst = {"Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",\
                                                    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", \
                                                    'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status',\
                                                    'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID',\
                                                    'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position'	,\
                                                    'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE','CANONICAL',\
                                                    'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', "SAS_MAF",\
                                                    "AA_MAF", 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO',\
                                                    'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'tcn', 'lcn', 'cf', 'purity', \
                                                    't_var_freq', 'expected_alt_copies', 'ccf_Mcopies', 'ccf_Mcopies_lower', 'ccf_Mcopies_upper', 'ccf_Mcopies_prob95', 'ccf_Mcopies_prob90', 'ccf_1copy', 'ccf_1copy_lower', 'ccf_1copy_upper', \
                                                    'ccf_1copy_prob95', 'ccf_1copy_prob90', 'ccf_expected_copies', 'ccf_expected_copies_lower', 'ccf_expected_copies_upper', 'ccf_expected_copies_prob95', 'ccf_expected_copies_prob90', 'clonality',\
                                                    'facets_fit', 'reviewer_set_purity', 'use_only_purity_run', 'use_edited_cncf', 'cncf_file_used'}
                                        if header== headtst:
                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 


                                elif fi.endswith(".ccf.maf"):
                                    # ccf check
                                    must_have_files.remove(".ccf.maf")
                                    with open(samplepath+"/"+runfi+"/"+fi) as segfile:
                                        for line in segfile:
                                            header=set(line.replace("\n","").split("\t"))
                                            break
                                        headtst = {'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1',	\
                                                    'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2',\
                                                    'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', \
                                                    'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', \
                                                    'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP',\
                                                    'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF',\
                                                    'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO',\
                                                    'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'flanking_bps', 'variant_id', 'variant_qual', 'ExAC_AF_Adj', \
                                                    'ExAC_AC_AN_Adj', 'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE', 'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER', 'Caller', 'COMMENTS', 'mutation_effect',\
                                                    'oncogenic', 'LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'LEVEL_R3', 'Highest_level', 'citations', 'oncokb_version', 'oncokb_version_date', 'residue', 'start_residue', 'end_residue',\
                                                    'variant_length', 'snv_hotspot', 'threeD_hotspot', 'indel_hotspot_type', 'indel_hotspot', 'Hotspot', 'driver', 'tcn', 'lcn', 'cf', 'purity', 't_var_freq', 'expected_alt_copies', 'ccf_Mcopies', 'ccf_Mcopies_lower',\
                                                    'ccf_Mcopies_upper', 'ccf_Mcopies_prob95', 'ccf_Mcopies_prob90', 'ccf_1copy', 'ccf_1copy_lower', 'ccf_1copy_upper', 'ccf_1copy_prob95', 'ccf_1copy_prob90', 'ccf_expected_copies', 'ccf_expected_copies_lower', \
                                                    'ccf_expected_copies_upper', 'ccf_expected_copies_prob95', 'ccf_expected_copies_prob90', 'clonality', 'facets_fit', 'reviewer_set_purity', 'use_only_purity_run', 'use_edited_cncf', 'cncf_file_used'}
                                        if header== headtst:
                                            pass
                                        else:
                                            if len(header)> len(headtst):
                                                print(samplepath+"/"+runfi+"/"+fi)
                                                print(header.difference(headtst))
                                                print("ccf.maf")
                                            patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi+"/"+fi) #Files missing cols 

                        else:
                            
                            if len(must_have_files) !=0:
                                for m in must_have_files:
                                    #Missing files 
                                    if m != "purity.seg" or m !=  "hisens.seg":
                                        patientIDorganized[samplefolder][2].append(samplepath+runfi+"/*"+m)
                                    
                        

                    else:


                        if os.path.getsize(samplepath+"/"+runfi) == 0 :
                            #Check file size and append to empty files
                                patientIDorganized[samplefolder][0].append(samplepath+"/"+runfi)

                        if runfi == "facets_qc.txt":
                            must_have_outerfiles.remove("facets_qc.txt")
                            with open(samplepath+"/"+runfi) as segfile:
                                for line in segfile:
                                    header=set(line.replace("\n","").split("\t"))
                                    break
                                headtst = {'tumor_sample_id', 'path', 'fit_name', 'purity_run_version', 'purity_run_prefix', 'purity_run_Seed', 'purity_run_cval', 'purity_run_nhet', 'purity_run_snp_nbhd', \
                                            'purity_run_ndepth', 'purity_run_Purity', 'purity_run_Ploidy', 'purity_run_dipLogR', 'purity_run_alBalLogR', 'hisens_run_version', 'hisens_run_prefix', 'hisens_run_Seed', \
                                            'hisens_run_cval', 'hisens_run_nhet', 'hisens_run_snp_nbhd', 'hisens_run_ndepth', 'hisens_run_hisens', 'hisens_run_Purity', 'hisens_run_Ploidy', 'hisens_run_dipLogR',\
                                            'manual_note', 'is_best_fit', 'purity', 'ploidy', 'dipLogR', 'dipLogR_flag', 'n_alternative_dipLogR', 'wgd', 'fga', 'n_dip_bal_segs', 'frac_dip_bal_segs', 'n_dip_imbal_segs', \
                                            'frac_dip_imbal_segs', 'n_amps', 'n_homdels', 'frac_homdels', 'n_homdels_clonal', 'frac_homdels_clonal', 'n_cn_states', 'n_segs', 'n_cnlr_clusters', 'n_lcn_na', 'n_loh', \
                                            'frac_loh', 'n_segs_subclonal', 'frac_segs_subclonal', 'n_segs_below_dipLogR', 'frac_below_dipLogR', 'n_segs_balanced_odd_tcn', 'frac_balanced_odd_tcn', 'n_segs_imbalanced_diploid_cn', \
                                            'frac_imbalanced_diploid_cn', 'n_segs_lcn_greater_mcn', 'frac_lcn_greater_mcn', 'n_snps', 'n_het_snps', 'frac_het_snps', 'n_snps_with_300x_in_tumor', 'n_het_snps_with_300x_in_tumor', \
                                            'n_het_snps_hom_in_tumor_1pct', 'n_het_snps_hom_in_tumor_5pct', 'frac_het_snps_hom_in_tumor_1pct', 'frac_het_snps_hom_in_tumor_5pct', 'mean_cnlr_residual', 'sd_cnlr_residual',\
                                            'n_segs_discordant_tcn', 'frac_discordant_tcn', 'n_segs_discordant_lcn', 'frac_discordant_lcn', 'n_segs_discordant_both', 'frac_discordant_both', 'n_segs_icn_cnlor_discordant', \
                                            'frac_icn_cnlor_discordant', 'mafr_median_all', 'mafr_median_clonal', 'mafr_n_gt_1', 'facets_suite_version', 'facets_qc_version', 'homdel_filter_pass', 'diploid_bal_seg_filter_pass', \
                                            'diploid_imbal_seg_filter_pass', 'waterfall_filter_pass', 'hyper_seg_filter_pass', 'high_ploidy_filter_pass', 'valid_purity_filter_pass', 'diploid_seg_filter_pass', \
                                            'em_cncf_icn_discord_filter_pass', 'dipLogR_too_low_filter_pass', 'subclonal_genome_filter_pass', 'icn_allelic_state_concordance_filter_pass', 'contamination_filter_pass', 'facets_qc'}
                                if header== headtst:
                                    pass
                                else:
                                    if len(header)> len(headtst):
                                        print(samplepath+"/"+runfi)
                                        print(header.difference(headtst))
                                        print("facets_qc.txt")
                                    patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi) #Files missing cols 
                        if runfi.endswith("nonsignedout.maf"):
                            must_have_outerfiles.remove("nonsignedout.maf")
                            with open(samplepath+"/"+runfi) as segfile:
                                for line in segfile:
                                    header=set(line.replace("\n","").split("\t"))
                                    break
                                headtst = {'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER'}

                                if header== headtst:
                                    pass

                                else:
                                    if len(header)> len(headtst):
                                        print(samplepath+"/"+runfi)
                                        print(header.difference(headtst))
                                        print("nonsignedout.maf")
                                    patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi) #Files missing cols 
                        
                        elif runfi.endswith("maf"):
                            must_have_outerfiles.remove(".maf")
                            with open(samplepath+"/"+runfi) as segfile:
                                for line in segfile:
                                    header=set(line.replace("\n","").split("\t"))
                                    break
                                headtst = {'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO', 'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'flanking_bps', 'variant_id', 'variant_qual', 'ExAC_AF_Adj', 'ExAC_AC_AN_Adj', 'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE', 'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER', 'Caller', 'COMMENTS', 'mutation_effect', 'oncogenic', 'LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'LEVEL_R3', 'Highest_level', 'citations', 'oncokb_version', 'oncokb_version_date', 'residue', 'start_residue', 'end_residue', 'variant_length', 'snv_hotspot', 'threeD_hotspot', 'indel_hotspot_type', 'indel_hotspot', 'Hotspot', 'driver'}
                                if header== headtst: 

                                    pass

                                else:
                                    if len(header)> len(headtst):
                                        print(samplepath+"/"+runfi)
                                    patientIDorganized[samplefolder][3].append(samplepath+"/"+runfi)  #Files missing cols 

                        #checks files in sample folder
            except:
                print("Error, permissions?")
                continue        

            #     print("not")

            patientIDorganized[samplefolder][1] = must_have_outerfiles #Missing files 
            
            patientIDorganized[samplefolder][6]= fit_name_li      #Missing fit folders
        # except:
        #     # print (bcolors.FAIL)
        #     print ("\t\tError in integrity_check. Terminating execution.")
        #     continue
        # print (e)
        # print (bcolors.ENDC)



    if EvalDirectory:

        outfile="ID\temptyfiles\tmissing_filesin_samplefolder\tMissingfilesinFitfolder\tmissingcols\ttumorbam_size_pass\tnormalbam_size_pass\tmanifest_file_labelerror\tmanifest_folder_missing\tConsent\tFitCount\tmultiple_in_manifest\tIn_Clinical_file\n"
        for keyID in patientIDorganized: 
            outfile+= keyID + "\t" + ",".join(patientIDorganized[keyID][0])+ "\t" + ",".join(patientIDorganized[keyID][1])+ "\t" + ",".join(patientIDorganized[keyID][2]) + "\t" + ",".join(patientIDorganized[keyID][3]) + "\t" + str(patientIDorganized[keyID][4][1]) + '\t' + str(patientIDorganized[keyID][4][0]) + "\t" + ",".join(patientIDorganized[keyID][5]) + "\t" + ",".join(patientIDorganized[keyID][6]) + "\t"+str(patientIDorganized[keyID][7]) + "\t"+str(patientIDorganized[keyID][8]) + "\t" + str(patientIDorganized[keyID][9]) + "\t" + str(patientIDorganized[keyID][10]) + "\n"

    else:
        outfile="ID\temptyfiles\tmissing_filesin_samplefolder\tMissingfilesinFitfolder\tmissingcols\ttumorbam_size_pass\tnormalbam_size_pass\tmanifest_file_labelerror\tmanifest_folder_missing\tConsent\tFitCount\tmultiple_in_manifest\n"
        for keyID in patientIDorganized: 
            outfile+= keyID + "\t" + ",".join(patientIDorganized[keyID][0])+ "\t" + ",".join(patientIDorganized[keyID][1])+ "\t" + ",".join(patientIDorganized[keyID][2]) + "\t" + ",".join(patientIDorganized[keyID][3]) + "\t" + str(patientIDorganized[keyID][4][1]) + '\t' + str(patientIDorganized[keyID][4][0]) + "\t" + ",".join(patientIDorganized[keyID][5]) + "\t" + ",".join(patientIDorganized[keyID][6]) + "\t"+str(patientIDorganized[keyID][7]) + "\t"+str(patientIDorganized[keyID][8])+ "\t" + str(patientIDorganized[keyID][9])+ "\n"


    date = datetime.now()
    datestr = str(date.isoformat(sep = "_",timespec='minutes')).replace(":","-")
    
    with open("FACETSintegrity_{datest}.tsv".format(datest= datestr),'w') as outf:
        outf.writelines(outfile)

    print("Sample number:" + str(sample_num))
    print("Fit number:" + str(fit_num))
    print("Fit Stats (reviewed_best_fit):" + str(fitstats[0]))
    print("Fit Stats (acceptable):" + str(fitstats[1]))
    print("Fit Stats (default):" + str(fitstats[2]))
    print("Fit Stats (not_reviewed):" + str(fitstats[3]))
    print("Fit Stats (No fit):" + str(fitstats[4]))
    print("Many normals vs same tumor number:" + str(manyNormal_sameTumor))


def test_facetsMeta(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    # facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    # prepared_metadata.buildFacetsMeta()
    # test_dataset = FacetsDataset(prepared_metadata)
    # #test_dataset.setCancerTypeFilter(["Lung"])
    # test_dataset.buildFacetsDataset()
    # integrity_check(test_dataset, EvalDirectory = False)
    # test_dataset.setPurityFilter(.5,1)
    integrity_check(inFacetsMeta = prepared_metadata)
    # integrity_check(inFacetsDataset = test_dataset)








if __name__ == '__main__':
    test_facetsMeta(True, True)

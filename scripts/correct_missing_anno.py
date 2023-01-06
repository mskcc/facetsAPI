import sys
import os

#change this to wherever the facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')
# sys.path.insert(1, '/juno/work/ccs/orgeraj/facetsAPI')

from facetsAPI import *

#This function will iterate a facets directory structure and check that mafs exist in each fit directory.
#If a maf is missing anywhere, then the maf generation and facetsPreview::generate_genomic_annotations() function
#will be run for the sample.
def correct_missing_annotations(useSingleRun, allowDefaults):
    clinical_sample_file  = ""#"/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    #facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    
    prepared_metadata.selectSamplesFromFile("/juno/work/ccs/pricea2/pipelines/facetsAPI/tests/test_samples.txt")

    prepared_metadata.buildFacetsMeta()

    fp_config = "/juno/work/ccs/bandlamc/git/ccs-cron/impact_facets/config_facets_preview.json"
    cbio_maf = "/work/ccs/shared/resources/impact/cbio_mutations/data_mutations_extended.oncokb.vep.maf"
    cbio_nonsigned_maf = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/cbio_mutations/data_mutations_nonsignedout.vep.maf"
    
    fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf)
    fp_tools.loadModule("R/R-3.6.3")

    #Sometimes you will get multiple fits, that will cause the job to queue multiple times.
    #Because facetspreview::generate_genomic_annotations will run at the sample base level on all fits
    #we don't want to queue multiple times for the same folder. So get unique paths before submitting.
    unique_paths = {} 
    for key in prepared_metadata.master_file_dict:
        cur_files = prepared_metadata.master_file_dict.get(key)
        cur_long_id = key.split('#')[0]
        cur_short_id = cur_long_id.split('_')[0]
        cur_ccf_maf = cur_files[7]
        cur_nonsigned_maf = cur_files[8]
        cur_path = os.path.dirname(cur_files[3]) + "/"
        if cur_path not in unique_paths:
            if cur_ccf_maf == "" or cur_nonsigned_maf == "":
                unique_paths[cur_path] = [cur_short_id, cur_long_id, cur_path, cur_ccf_maf]

    print(str(len(prepared_metadata.master_file_dict)) + " items in master file dict.")
    print(str(len(unique_paths)) + " unique paths need correction.")

    for cur_unique_path in unique_paths:
        cur_correct = unique_paths.get(cur_unique_path)
        print("\t\t\tRunning FacetsPreview generate_genomic_annotations() for " + cur_correct[1])
        fp_tools.runGenerateGenomicAnnotations(cur_correct[0], cur_correct[1], cur_correct[2])

    print("\t\tCompleted correct_missing_annotations.")

if __name__ == '__main__':
    correct_missing_annotations(False, False)
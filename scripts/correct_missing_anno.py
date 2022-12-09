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
    prepared_metadata.buildFacetsMeta()

    fp_config = "/juno/work/ccs/bandlamc/git/ccs-cron/impact_facets/config_facets_preview.json"
    cbio_maf = "/work/ccs/shared/resources/impact/cbio_mutations/data_mutations_extended.oncokb.vep.maf"
    cbio_nonsigned_maf = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/cbio_mutations/data_mutations_nonsignedout.vep.maf"
    
    fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf)
    fp_tools.loadModule("R/R-3.6.3")

    for key in prepared_metadata.master_file_dict:
        cur_files = prepared_metadata.master_file_dict.get(key)
        cur_long_id = key.split('#')[0]
        cur_short_id = cur_long_id.split('_')[0]
        cur_ccf_maf = cur_files[7]
        cur_nonsigned_maf = cur_files[8]
        cur_path = os.path.dirname(cur_files[3]) + "/"
    
        if cur_ccf_maf == "" or cur_nonsigned_maf == "":
            print("\t\t\tRunning FacetsPreview generate_genomic_annotations() for " + key)
            fp_tools.runGenerateGenomicAnnotations(cur_short_id, cur_long_id, cur_path)

    print("\t\tCompleted correct_missing_annotations.")

if __name__ == '__main__':
    correct_missing_annotations(False, False)
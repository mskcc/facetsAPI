import sys
import os
import shutil

#change this to wherever the facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/juno/work/ccs/shared/resources/impact/knowledge-systems/2024-10-29/dmp_pipeline_annotated/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity", "/juno/work/ccs/pricea2/tmp/neoadjuvant/round2/facetsAPI.meta")
    
    #We just want to look at a single run per sample, looking for best fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)
    
    #Read in the list of IDs we are selecting from a file.
    prepared_metadata.selectSamplesFromFile("/juno/work/ccs/pricea2/tmp/neoadjuvant/round2/partA_consented.txt")
    
    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    print(prepared_metadata.fit_map)

    #Build our FacetsDataset Object.
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    #test_dataset.copyDatasetToFolder("/juno/work/ccs/pricea2/tmp/neoadjuvant/")

    ext_tools = ExtTools(prepared_metadata)
    #ext_tools.makeMergedFile(MetaDictMap.UNADJUSTED_SEG_FILE, "/juno/work/ccs/pricea2/pipelines/facetsAPI/tests/ascet_test/merged_seg.txt")
    ext_tools.runAscets("/juno/work/ccs/pricea2/tmp/neoadjuvant/round2/", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=False)
    test_dataset.writeReport("/juno/work/ccs/pricea2/tmp/neoadjuvant/round2/neoadjuvant_sample_report.txt")

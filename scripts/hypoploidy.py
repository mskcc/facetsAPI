import sys
import os

#change this to wherever the facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/oncokb-annotated-msk-impact/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(True,False)

    prepared_metadata.selectSamplesFromFile("/juno/work/ccs/pricea2/pipelines/facetsAPI/scripts/test_hypo.txt")
    
    prepared_metadata.buildFacetsMeta()


    hypo_dataset = FacetsDataset(prepared_metadata)

    filter_on_qc_fail       = True     # If this is True, facets_qc failed samples will be excluded from analysis.
    purity_filter_thresh    = 0.2      # This is the level of purity at which we disqualify samples from consideration.
    hypo_dataset.setPurityFilter(purity_filter_thresh)
    #hypo_dataset.setFacetsQCFilter(filter_on_qc_fail)
    hypo_dataset.buildFacetsDataset()
    hypo_dataset.runAlterationAnalysis()

    #Write output.
    hypo_dataset.writeAlterationData("bestOnly_02purityFilter_hypoploidy_data.txt")
    hypo_dataset.writeReport("bestOnly_02purityFilter_hypoploidy_report.txt")

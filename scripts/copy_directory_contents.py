import sys
import os

# facetsAPI must be importable: `pip install -e /path/to/facetsAPI` (editable) or add the repository
# root to PYTHONPATH. (Earlier versions hardcoded a cluster path here.)

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = ""#"/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/oncokb-annotated-msk-impact/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/tcga/all3/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    
    #We just want to look at a single run per sample, looking for best fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)
    
    #Read in the list of IDs we are selecting from a file.
    #prepared_metadata.selectSamplesFromFile("/path/to/facetsAPI/scripts/DMPID_CSFonly.txt")
    
    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    #Build our FacetsDataset Object.
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    #Copy our selected samples to the new output location.
    test_dataset.writeReport("/path/to/facetsAPI/scripts/TCGA_report.txt")
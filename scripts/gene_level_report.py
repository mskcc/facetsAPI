import sys
import os

# Use the facetsAPI in this checkout (the folder above scripts/), not any installed copy.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/facets/all/"
    input_samples         = "/path/to/input_samples.txt"   # one tumor sample id per line (P-0000000-T01-IM6)
    output_dir            = "/path/to/output_dir/"
    target_genes          = ["TP53"]

    os.makedirs(output_dir, exist_ok=True)

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")

    #We just want to look at a single run per sample, looking for best fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)

    #Read in the list of IDs we are selecting from a file.
    prepared_metadata.selectSamplesFromFile(input_samples)

    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    #Build our FacetsDataset Object.
    test_dataset = FacetsDataset(prepared_metadata)
    #test_dataset.setFacetsQCFilter(True)
    test_dataset.buildFacetsDataset()

    #Gene level copy number (tcn/lcn/cf/cn_state...) for our target genes. These come from gene_level.txt (hisens run).
    test_dataset.writeGeneLevelReport(os.path.join(output_dir, "gene_level_report.txt"), target_genes)

    #Write a report of the selected samples (fit, purity, ploidy, WGD, FGA, QC...).
    test_dataset.writeReport(os.path.join(output_dir, "report.txt"))

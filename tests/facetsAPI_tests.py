import sys

#change this to wherever the facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')
# sys.path.insert(1, '/juno/work/ccs/orgeraj/facetsAPI')

from facetsAPI import *

def test_alterationFunctions(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()
    test_dataset.runAlterationAnalysis()
    test_dataset.writeAlterationData("test_alterations.txt")

#Test building a facetsDataset object.  If useSingleRun is True, the best fit/acceptable fits will be selected.
#Otherwise, an object that considers all runs for each sample will be generated.  
def test_facetsDataset(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(facets_repo_path =facets_dir, hisens_vs_purity = "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()
    test_dataset = FacetsDataset(prepared_metadata)
    
    test_dataset.buildFacetsDataset()
    test_dataset.printFacetsSampleById("P-0082589-T01-IM7_P-0082589-N01-IM7")
    sample = test_dataset.sampleList.get("P-0082589-T01-IM7_P-0082589-N01-IM7")
    # print(sample)
    # test_run = sample.runs[0]
    # test_segment = test_run.segments[0]
    # test_segment.printSegment()
    # test_genes = test_run.genes[0]
    # test_genes.printGene()
    

#Test building a FacetsMeta object.  If useSingleRun is True, the best fit/acceptable fits will be selected.
#Otherwise, an object that considers all runs for each sample will be generated.  
def test_facetsMeta(useSingleRun, allowDefaults):
    clinical_sample_file  = ""#"/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    #facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()
    test_dataset = FacetsDataset(prepared_metadata)
    #test_dataset.setCancerTypeFilter(["Lung"])
    test_dataset.buildFacetsDataset()

def test_facetsPurityByCF():
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "hisens")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.buildFacetsMeta()

    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    test_dataset.writePurityCFs("cfPurity_hisens_base.txt", True)
    test_dataset.writePurityCFs("cfPurity_hisens_cfEm.txt", False)


if __name__ == '__main__':

    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
    test_facetsMeta(False, False)
    #test_facetsDataset(True, True)
    #test_alterationFunctions(False, False)
    #test_facetsPurityByCF()
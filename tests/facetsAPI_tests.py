import sys
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')

from facetsAPI import *

def test_alterationFunctions(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()
    test_dataset = FacetsDataset()
    test_dataset.buildFacetsDataset(prepared_metadata)
    test_dataset.runAlterationAnalysis()
    test_dataset.writeAlterationData("test_alterations.txt")

#Test building a facetsDataset object.  If useSingleRun is True, the best fit/acceptable fits will be selected.
#Otherwise, an object that considers all runs for each sample will be generated.  
def test_facetsDataset(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()
    test_dataset = FacetsDataset()
    test_dataset.buildFacetsDataset(prepared_metadata)
    test_dataset.printFacetsSampleById("P-0000004-T01-IM3_P-0000004-N01-IM3")
    

#Test building a FacetsMeta object.  If useSingleRun is True, the best fit/acceptable fits will be selected.
#Otherwise, an object that considers all runs for each sample will be generated.  
def test_facetsMeta(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()


def test_facetsPurityByCF():
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "hisens")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.buildFacetsMeta()

    test_dataset = FacetsDataset()
    test_dataset.buildFacetsDataset(prepared_metadata)

    test_dataset.writePurityCFs("cfPurity_hisens_base.txt", True)
    test_dataset.writePurityCFs("cfPurity_hisens_cfEm.txt", False)


if __name__ == '__main__':

    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #test_facetsMeta(False, False)
    #test_facetsDataset(False, False)
    test_alterationFunctions(False, False)
    #test_facetsPurityByCF()
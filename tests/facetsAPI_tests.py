import sys
import os

# facetsAPI must be importable: `pip install -e /path/to/facetsAPI` (editable) or add the repository
# root to PYTHONPATH. (Earlier versions hardcoded a cluster path here.)

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
    #facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.selectSamplesFromFile("/path/to/facetsAPI/tests/test_samples.txt")
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

def test_tools(useSingleRun, allowDefaults):
    clinical_sample_file  = ""#"/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    #facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()

    fp_config = "/path/to/config_facets_preview.json"
    cbio_maf = "/work/ccs/shared/resources/impact/cbio_mutations/data_mutations_extended.oncokb.vep.maf"
    cbio_nonsigned_maf = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/cbio_mutations/data_mutations_nonsignedout.vep.maf"
    
    fp_tools = Tools(fp_config, cbio_maf, cbio_nonsigned_maf)
    # R is reached through FacetsConfig: set FACETSAPI_R_PREFIX (e.g. the IRIS spack module-use + module load
    # line) and FACETSAPI_R_LIB_LOC, or pass config=FacetsConfig(...) to FPTools. See doc/configuration.md.

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


        #print(cur_files)
        #print(cur_long_id)
        #print(cur_files[7])
        #print(cur_files[8])
        #sys.exit()

def test_ascets(useSingleRun, allowDefaults):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    #facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()

    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()
    
    ext_tools = ExtTools(prepared_metadata)
    #ext_tools.makeMergedFile(MetaDictMap.UNADJUSTED_SEG_FILE, "/path/to/facetsAPI/tests/ascet_test/merged_seg.txt")
    #ext_tools.runAscets("/path/to/facetsAPI/tests/ascet_test2/")
    test_dataset.writeReport("/path/to/facetsAPI/tests/ascet_test2/sample_report.txt")

if __name__ == '__main__':
    #clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    #facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
    #test_tools(False, False)
    #test_facetsMeta(False, False)
    test_ascets(True,True)
    #test_facetsDataset(True, True)
    #test_alterationFunctions(False, False)
    #test_facetsPurityByCF()
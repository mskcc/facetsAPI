import sys
import os

#change this to wherever the facetsAPI is stored
sys.path.insert(1, '/juno/work/ccs/pricea2/pipelines/facetsAPI')
# sys.path.insert(1, '/juno/work/ccs/orgeraj/facetsAPI')

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/oncokb-annotated-msk-impact/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    
    #We just want to look at a single run per sample, looking for best fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)
    
    #Read in the list of IDs we are selecting from a file.
    prepared_metadata.selectSamplesFromFile("/juno/work/ccs/pricea2/pipelines/facetsAPI/tests/test_samples.txt")
    
    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    #Build our FacetsDataset Object.
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    #Now we have all the info we want in one place, so write a custom function that outputs the report we want.
    #Write the info we want out to a file.  Keep track of missing samples so we can report as well.
    with open("sam_cohort_report2.txt", 'w') as outfile:
        outfile.write("ID\tfitDir\tfitType\tCancer Type\tCancer Type Detail\tFacets Purity\tClinical Purity\tOnkoCode\tPloidy\tdipLogR\tcVal\tWGD\tFGA\tfrac_loh\tfacets_qc_pass\ttmb\tmsi" + "\n")
        found_list = []
        for cur_run in test_dataset.runList:
            print(cur_run.id)
            cur_tumor_id = cur_run.id[0:17]
            cur_fit_item = prepared_metadata.fit_map.get(cur_run.id)
            print("fit " + cur_fit_item[0])
            print("wgd: " + str(cur_run.wgd))

            outfile.write(str(cur_run.id) + "\t")
            outfile.write(str(cur_run.fitDir) + "\t")
            outfile.write(str(cur_fit_item[0]) + "\t")
            outfile.write(str(cur_run.cancerType) + "\t")
            outfile.write(str(cur_run.cancerTypeDetail) + "\t")
            outfile.write(str(cur_run.purity) + "\t")
            outfile.write(str(cur_run.clinicalPurity) + "\t")
            outfile.write(str(cur_run.onkoCode) + "\t")
            outfile.write(str(cur_run.ploidy) + "\t")
            outfile.write(str(cur_run.dipLogR) + "\t")
            outfile.write(str(cur_run.cval) + "\t")
            outfile.write(str(cur_run.wgd) + "\t")
            outfile.write(str(cur_run.fga) + "\t")
            outfile.write(str(cur_run.frac_loh) + "\t")
            outfile.write(str(cur_run.facets_qc) + "\t")
            outfile.write(str(cur_run.tmb) + "\t")
            outfile.write(str(cur_run.msi) + "\t")
            outfile.write("\n")

            if cur_tumor_id in prepared_metadata.samples_from_file:
                #print(cur_tumor_id + " found okay.")
                found_list.append(cur_tumor_id)

        for expected in prepared_metadata.samples_from_file:
            if expected not in found_list:
                outfile.write(expected + "\t" + "MISSING\n")

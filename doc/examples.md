## Examples

#### Selecting a FACETS data subset.

The following code will build a FacetsDataset with a set of filters applied so that only samples passing FACETS QC, 
who have been specified as "best fits", and with a minimum purity cutoff of 0.2.  After selecting and building this dataset the user could
manipulate or output it as they require.

```
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, True, "purity", "all_meta.dat")

    hypo_dataset = FacetsDataset()
    filter_on_qc_fail       = True     # If this is True, facets_qc failed samples will be excluded from analysis.
    purity_filter_thresh    = 0.2      # This is the level of purity at which we disqualify samples from consideration.
    hypo_dataset.setPurityFilter(purity_filter_thresh)
    hypo_dataset.setFacetsQCFilter(filter_on_qc_fail)
    hypo_dataset.buildFacetsDataset(prepared_metadata)

```

#### Calculating aneuploidy scores for FACETS data using Ascets via FacetsAPI.

[Ascets](https://github.com/beroukhim-lab/ascets) (Arm-level Somatic Copy-number Events in Targeted Sequencing) is a tool for calculating aneuploidy scores for arm-level data.  Several MSK researchers have been interested in seeing the results of ascets based on FACETS data.  The following code is an example of using FacetsAPI to execute such a run:

```
import sys
import os

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/oncokb-annotated-msk-impact/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    
    #We just want to look at a single run per sample, looking for best fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)
    
    #Read in the list of IDs we are selecting from a file.
    prepared_metadata.selectSamplesFromFile("/path/to/input_sample_list.txt")
    
    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    #Build our FacetsDataset Object.
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    ext_tools = ExtTools(prepared_metadata)
    ext_tools.runAscets("/path/to/output/", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True)
    test_dataset.writeReport("/path/to/output/sample_report.txt")

```

#### Performing hypodiploidy identification analysis.

Work to identify and investigate patterns of LOH in FACETS data is being performed at MSK.  The following code is an example of using FacetsAPI to generate a curated list of samples and sample-level calls for hypodiploidy samples.  This code splits FACETS segment data by chromosome arm and identifies all samples with a hypodiploidy as defined by provided threshold parameters.  See [FacetsDataset](facetsdataset.md)'s runAlterationAnalysis() for additional details.

```
import sys
import os

from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/knowledge-systems/11-14-22/oncokb-annotated-msk-impact/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(True,True)
    prepared_metadata.buildFacetsMeta()

    hypo_dataset = FacetsDataset(prepared_metadata)
    filter_on_qc_fail       = True     # If this is True, facets_qc failed samples will be excluded from analysis.
    purity_filter_thresh    = 0.2      # This is the level of purity at which we disqualify samples from consideration.
    hypo_dataset.setPurityFilter(purity_filter_thresh)
    hypo_dataset.setFacetsQCFilter(filter_on_qc_fail)
    hypo_dataset.buildFacetsDataset()
    hypo_dataset.runAlterationAnalysis()

    #Write output.
    hypo_dataset.writeAlterationData("all_hypoploidy.txt")

```

### Produce a report with information from a select set of samples.

A researcher wanted to look at several pieces of data for a specific set of FACETS samples.  The following code accepts a sample list from a file and compiles an output report with the relevant information for the requested samples.

```
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
    
    #We just want to look at a single run per sample, looking for best/acceptable fits. Default is acceptable if not.
    prepared_metadata.setSingleRunPerSample(True,True)
    
    #Read in the list of IDs we are selecting from a file. One sample per line.
    prepared_metadata.selectSamplesFromFile("/path/to/input_test_samples.txt")
    
    #Build our FacetsMeta Object.
    prepared_metadata.buildFacetsMeta()

    #Build our FacetsDataset Object and write a report to file.
    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()
    test_dataset.writeReport("/path/to/report.txt")
```


### Produce genomic annotation for Facets Samples that are missing genomic annotation files.

At one point, it was found that some FACETS samples in the impact repository were missing genomic annotation.  The following code is an example of using FacetAPI to execute FacetsPreview::generate_genomic_annotation() on the affected samples through the Juno LSF cluster.

```
import sys
import os

from facetsAPI import *

#This function will iterate a facets directory structure and check that mafs exist in each fit directory.
#If a maf is missing anywhere, then the maf generation and facetsPreview::generate_genomic_annotations() function
#will be run for the sample.
def correct_missing_annotations(useSingleRun, allowDefaults):
    clinical_sample_file  = ""
    facets_dir            = "/my/facets/folder/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    
    prepared_metadata.buildFacetsMeta()

    fp_config = "/my/facets/config_facets_preview.json"
    cbio_maf = "/my/facets/data_mutations_extended.oncokb.vep.maf"
    cbio_nonsigned_maf = "/my/facets/data_mutations_nonsignedout.vep.maf"
    
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
```

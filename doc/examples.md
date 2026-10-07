## Examples

#### Selecting a FACETS data subset.

The following code will build a FacetsDataset with a set of filters applied so that only samples passing FACETS QC,
with one fit per sample (human-reviewed fits first, then autoQC fits, then the default fit), and with a minimum purity cutoff of 0.2.
After selecting and building this dataset the user could manipulate or output it as they require.

```python
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward and persist it to all_meta.dat.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity", "all_meta.dat")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.buildFacetsMeta()

    hypo_dataset = FacetsDataset(prepared_metadata)
    filter_on_qc_fail       = True     # If this is True, facets_qc failed samples will be excluded from analysis.
    purity_filter_thresh    = 0.2      # This is the level of purity at which we disqualify samples from consideration.
    hypo_dataset.setPurityFilter(purity_filter_thresh)
    hypo_dataset.setFacetsQCFilter(filter_on_qc_fail)
    hypo_dataset.buildFacetsDataset()
```

#### Choosing which fits count: reviewed only, autoQC only, default only.

```python
from facetsAPI import FacetsMeta, FacetsDataset, FitPolicy

meta = FacetsMeta(clinical_sample_file, facets_dir, "purity")
meta.setFitPolicy(FitPolicy.reviewed_only())    # or FitPolicy.autoqc_only(), FitPolicy.default_only(),
meta.setSingleRunPerSample(True, False)         #    FitPolicy(("autoqc", "reviewed", "default")) ...
meta.buildFacetsMeta()
print(meta.build_counts)                        # how many Best / Acceptable / AutoQC_Best / AutoQC_Pass / Default fits were chosen
```

#### Working with a 2n (facets2n) repository.

2n pairs carry a `clinical/` and a `research/` fit set. Pick one with `fit_class`, or load both as two runs of the same sample.

```python
from facetsAPI import *

meta = FacetsMeta(clinical_sample_file, "/path/to/impact_2n/facets/all/", "purity", fit_class="both")
meta.setSingleRunPerSample(True, True)
meta.buildFacetsMeta()
ds = FacetsDataset(meta)
ds.buildFacetsDataset()
for sample_id, sample in ds.sampleList.items():            # keyed by <tumor>_<normal>
    for run in sample.runs:                                 # one clinical + one research run
        print(sample_id, run.fit_class, run.fit_name, run.review_status, run.purity)
```

#### Calculating aneuploidy scores for FACETS data using Ascets via FacetsAPI.

[Ascets](https://github.com/beroukhim-lab/ascets) (Arm-level Somatic Copy-number Events in Targeted Sequencing) is a tool for calculating aneuploidy scores for arm-level data.  Several MSK researchers have been interested in seeing the results of ascets based on FACETS data.  The following code is an example of using FacetsAPI to execute such a run. R is reached through [FacetsConfig](configuration.md); on IRIS set `FACETSAPI_R_PREFIX` to the spack `module use ... && module load r/3.6.3` line (or pass `config=`). The run executes inline; submit the script itself with `sbatch` if you want it on a compute node.

```python
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.selectSamplesFromFile("/path/to/input_sample_list.txt")
    prepared_metadata.buildFacetsMeta()

    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()

    r_config  = FacetsConfig(r_prefix="module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3")
    ext_tools = ExtTools(prepared_metadata, config=r_config)
    ext_tools.runAscets("/path/to/output/", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True)
    test_dataset.writeReport("/path/to/output/sample_report.txt")
```

#### Performing hypodiploidy identification analysis.

Work to identify and investigate patterns of LOH in FACETS data is being performed at MSK.  The following code is an example of using FacetsAPI to generate a curated list of samples and sample-level calls for hypodiploidy samples.  This code splits FACETS segment data by chromosome arm and identifies all samples with a hypodiploidy as defined by provided threshold parameters.  See [FacetsDataset](facetsdataset.md)'s runAlterationAnalysis() for additional details.

```python
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.buildFacetsMeta()

    hypo_dataset = FacetsDataset(prepared_metadata)
    hypo_dataset.setPurityFilter(0.2)
    hypo_dataset.setFacetsQCFilter(True)
    hypo_dataset.buildFacetsDataset()
    hypo_dataset.runAlterationAnalysis()
    hypo_dataset.writeAlterationData("all_hypoploidy.txt")
```

#### Produce a report with information from a select set of samples.

```python
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
    facets_dir            = "/path/to/impact/facets/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(True, True)
    prepared_metadata.selectSamplesFromFile("/path/to/input_test_samples.txt")   # one sample id per line
    prepared_metadata.buildFacetsMeta()

    test_dataset = FacetsDataset(prepared_metadata)
    test_dataset.buildFacetsDataset()
    test_dataset.writeReport("/path/to/report.txt")
```

#### Produce genomic annotation for Facets Samples that are missing genomic annotation files.

At one point, it was found that some FACETS samples in the impact repository were missing genomic annotation.  The following code is an example of using FacetsAPI to execute FacetsPreview::generate_genomic_annotation() on the affected samples. The R wrapper runs inline through [FacetsConfig](configuration.md).

```python
import os
from facetsAPI import *

def correct_missing_annotations(useSingleRun, allowDefaults):
    clinical_sample_file  = ""   # directory structure only; clinical data is not needed here
    facets_dir            = "/my/facets/folder/all/"

    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun, allowDefaults)
    prepared_metadata.buildFacetsMeta()

    fp_config = "/my/facets/config_facets_preview.json"
    cbio_maf = "/my/facets/data_mutations_extended.oncokb.vep.maf"
    cbio_nonsigned_maf = "/my/facets/data_mutations_nonsignedout.vep.maf"
    r_config = FacetsConfig(r_prefix="module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3",
                            r_lib_loc="/my/R_libs")
    fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf, config=r_config)

    # generate_genomic_annotations works at the sample (or 2n class) level, so run once per directory.
    unique_paths = {}
    for key in prepared_metadata.master_file_dict:
        cur_files = prepared_metadata.master_file_dict.get(key)
        cur_long_id = prepared_metadata.pairIdForKey(key)
        cur_short_id = FacetsMeta.shortId(key)
        cur_path = cur_files[MetaDictMap.SAMPLE_BASE_DIR] + "/"
        if cur_path not in unique_paths:
            if cur_files[MetaDictMap.CCF_MAF_FILE] == "" or cur_files[MetaDictMap.NONSIGNED_MAF_FILE] == "":
                unique_paths[cur_path] = [cur_short_id, cur_long_id, cur_path]

    print(str(len(unique_paths)) + " unique paths need correction.")
    for cur_correct in unique_paths.values():
        fp_tools.runGenerateGenomicAnnotations(cur_correct[0], cur_correct[1], cur_correct[2])

if __name__ == '__main__':
    correct_missing_annotations(False, False)
```

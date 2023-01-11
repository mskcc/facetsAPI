## FacetsDataset

### Initializing a FacetsDataset

After building a [FacetsMeta](facetsmeta.md) object, it is possible to define a FacetsDataset that represents all or a specified subset of facets data based on a variety of conditions.  A FacetsDataset contains a list of [FacetsSamples](facetssample.md) and/or [FacetsRuns](facetsrun.md).  When performing operations on a range of samples or runs, iterating the constructed lists contained in a FacetsDataset is the suggested approach.  The FacetsDataset object also contains functions that apply operations to the entire dataset, or that write output for a range of samples/runs.

The initialization function of a FacetsDataset class does not require any parameters. Instead, building a FacetsDataset object generally involves three steps.
* Create an empty FacetsDataset: `test_dataset = FacetsDataset()`
* Apply desired filter/selection functions.
* Execute `buildFacetsDataset() function.`

#### Selection / Filtration Functions
* setCancerTypeFilter(string[] selectedCancerTypes)
  * This function applies a filter to the FacetsDataset build so that only indicated cancer types will be selected. Input is expected to a list of strings, even if the length of the selection is 1. Note that cancer types are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setCancerTypeFilter(["Example Cancer Type 1", "Example Cancer Type 2"])`.
* setCancerTypeDetailFilter(string[] selectedCancerDetailTypes)
  * This function applies a filter to the FacetsDataset build so that only indicated detail cancer types will be selected. Input is expected to a list of strings, even if the length of the selection is 1. Note that cancer types are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setCancerTypeDetailFilter(["Example Cancer Type 1", "Example Cancer Type 2"])`.
* setPurityFilter(minPurity, maxPurity=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples with purity values in a given range are selected.  minPurity is required, however maxPurity, unless otherwise indicated, defaults to 1.
* setClinicalPurityFilter(minPurity, maxPurity=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples with clinical purity values in a given range are selected. Clinical purity is the purity value extracted from the data_clinical_sample.txt file, whereas purity, as used in setPurityFilter(), references FACETS calculated purity. minPurity is required, however maxPurity, unless otherwise indicated, defaults to 1.
* setOnkoCodeFilter(selectedOnkoCodes)
  *  This function applies a filter to the FacetsDataset build so that samples with specified onko codes are included. Input is expected to a list of strings, even if the length of the selection is 1. Note that onko codes are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setOnkoCodeFilter(["LUAD","LUSC"])`.
* setPloidyFilter(minPloidy, maxPloidy=1.0)
  *  This function applies a filter to the FacetsDataset build so that only samples within a given ploidy range are selected. minPloidy and maxPloidy are both required.
* setDipLogRFilter(minDipLogR, maxDipLogR)
  *  This function applies a filter to the FacetsDataset build so that only samples within a given dipLogR range are selected. minDipLogR and maxDipLogR are both required.
* setCvalFilter(minC, maxC)
  * This function applies a filter to the FacetsDataset build so that only samples within a given cValue range are selected. minC and maxC are both required.
* setWGDFilter(keepWGD)
  * This function applies a filter to the FacetsDataset build indicating whether to select only WGD samples, or only non-WGD samples.  If this is set to True, WGD samples will be selected.  If set to False, only non-WGD samples will be selected.
* setFGAFilter(minFGA, maxFGA)
  * This function applies a filter to the FacetsDataset build so that only samples within a given FGA range are selected.  minFGA is required, however maxFGA, unless otherwise indicated, defaults to 1.
* setFracLohFilter(minLoh, maxLoh)
  * This function applies a filter to the FacetsDataset build so that only samples within a given frac_loh range are selected.  minLoh is required, however maxLoh, unless otherwise indicated, defaults to 1.
* setFacetsQCFilter(keepQCPass)
  * This function applies a filter to the FacetsDataset build indicating whether to select only samples that pass FACETS QC, or only those who fail.  If this is set to True, only FACECTS QC PASS samples will be selected. If set to False, only FACECTS QC FAIL samples will be selected.
* setTMBFilter(minTMB, maxTMB)
  * This function applies a filter to the FacetsDataset build so that only samples within a given TMB range are selected. minTMB and maxTMB are both required. 
* setMSIFilter(minMSI, maxMSI)
  * This function applies a filter to the FacetsDataset build so that only samples within a given MSI range are selected. minMSI and maxMSI are both required. 


#### Dataset Functions
* buildFacetsDataset(facets_meta_object)
  * This function will populate a FacetsDataset object with samples and/or runs that meet the specified selection criteria indicated by using the above filtration functions.  This function accepts a fully constructed [FacetsMeta](facetsmeta.md) object as a parameter, and uses the directory mappings in the [FacetsMeta](facetsmeta.md) object to navigate and parse the FACETS file system for target files.  This function has two relvant variables in the object that can be referenced for custom processes `sampleList` and `runList`. sampleList is python dictionary of the form {sample_id -> FacetsSample} and runList is an array of type FacetsRun.  Either of these structures contain the entire set of data, however they are structured differently.  A [FacetsSample](facetssample.md) represents a sample that may have multiple runs, whereas a [FacetsRun](facetsrun.md) represents a single run.  References to both of these structrues are stored to allow for differential processing options and do not contribute to increased memory, as FacetsSample and FacetsRun are accessed by reference, and not by value.
* runAlterationAnalysis()
  * This function will run alteration analysis on a constructed FacetsDataset object.  This process will run several functions in each [FacetsRun](facetsrun.md).  Specifically, for each sample, arms will be defined for all [FacetsSegments](facetssegment.md), then LoH and Gains will be calculated at the arm level, CF levels will be calculated, and run level alteration calls will be made.  (i.e. "Hypoploidy")
* writeAlterationData(outfile_path)
  * This function will write alteration data for all samples in the FacetsDataset to the specified outfile_path.  This includes arm level gain/loss calls, CF levels, and run-level calls.
* writePurityCFs(outfile_path, use_base_cf=False, use_long_id=True)
  * This function will calculate and write calculated purity data for each sample.  Purity is defined as the greatest CF value that is not NA or 1. If use_base_cf is true, cf will be used rather than cf.em.  If use_long_id is True, the long_id including the paired normal information will be used when writing IDs. 
* printFacetsSampleById(id)
  * This function will print a specific sample in this facets dataset. id is expected to be a string. 
* writeDatasetSummary()
  * This function will output a text file containing information tat summarizes a facets dataset. This includes metrics like number of best/acceptable fits, number of failed samples that couldn't build, number of cancer types, etc.  The output's name contains the date/time
* createHistogram(variable, inbins=10)
  * This function creates a histogram from the requested variable.  It defaults to 10 automated bins for numerical variables but that can be changed. For strings like cancerType or boolean values like WGD it uses however many bins are needed.
  * Available variables are as follows: purity, clinicalPurity, ploidy, dipLogR, cval, fga, frac_loh, tmb, msi, cancerType, oncoCode, wgd, and review_status. 


#### Examples

```python
  # Build a FacetsDataset that contains only Breast Cancer type samples with purity between 0.2 and 0.8. 
  
   prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
   prepared_metadata.setSingleRunPerSample(True, True)
   prepared_metadata.buildFacetsMeta()
   test_dataset = FacetsDataset()
   test_dataset.setPurityFilter(0.2, 0.8)
   test_dataset.setCancerTypeFilter(["Breast Cancer"])
   test_dataset.buildFacetsDataset(prepared_metadata)
```

```python
  # Build a FacetsDataset that contains a few onkocode specific samples with purity between 0.3 and 1, 
  # and run alteration analysis. Output the result to file. 
  
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, False, "hisens")
  prepared_metadata.buildFacetsMeta()
  test_dataset = FacetsDataset()
  test_dataset.setPurityFilter(0.3)
  test_dataset.setOnkoCodeFilter(["SCLC","CSCLC"])
  test_dataset.printFacetsSampleById("P-0000004-T01-IM3_P-0000004-N01-IM3")
  test_dataset.runAlterationAnalysis()
  test_dataset.writeAlterationData("test_alterations_data.txt")
```

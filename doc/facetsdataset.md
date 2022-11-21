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
  *  This function applies a filter to the FacetsDataset build so that samples with specified onko codes are included. Input is expected to a list of strings, even if the length of the selection is 1. Note that onko codes are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setOnkoCodeFilter(["LUAD","LUSC")`.
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


* printLogo() 
  * A simple method for printing a FACETS API logo. 
  **This function is run automatically during buildFacetsMeta()**.
* parseClinicalSample()
  * This function will accept a clinical sample file and scan each sample's corresponding facets directory.
  the facets_review.manifest file will be read to identify
  the appropriate fit to use in analysis.  A dictionary mapping
  samples to their appropriate directories is built in master_file_dict.
  This function will also build dictionaries mapping sample id to
  cancer type and to patient id. 
  **This function is run automatically during buildFacetsMeta()**.
* buildFacetsMeta()
  * This function will populate details of the FacetsMeta object with data from the data_clinical_sample file.  This function also handles data persistence.  If this FacetsMeta object was constructed using the optional persistence parameter, this function will attempt to load data from the specified file if it exists.  If the specified file does not exist, then this function will parse the clinical data and store the specified file for future loading.
* setVerbose(bool)
  * This function accepts True or False as a parameter, and will activate verbose mode if set to True.  This will make warning messages visible during runtime and show more details of data processing for various processes.  By default, verbose mode is set to False.   


#### Examples

```python
  # This will build a FacetsMeta object that includes default fits, selects "purity" type data, 
  # and persists the FacetsMeta object into a file called all_meta.dat
  
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, True, "purity", "all_meta.dat")
  prepared_metadata.buildFacetsMeta()
```


```python
  # This will build a FacetsMeta object that includes only curated fits, 
  # selects "hisens" type data, and does not persist data.
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, False, "hisens")
  prepared_metadata.buildFacetsMeta()
```

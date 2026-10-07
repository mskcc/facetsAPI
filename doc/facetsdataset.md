## FacetsDataset

### Initializing a FacetsDataset

After building a [FacetsMeta](facetsmeta.md) object, it is possible to define a FacetsDataset that represents all or a specified subset of facets data based on a variety of conditions.  A FacetsDataset contains a list of [FacetsSamples](facetssample.md) and/or [FacetsRuns](facetsrun.md).  When performing operations on a range of samples or runs, iterating the constructed lists contained in a FacetsDataset is the suggested approach.  The FacetsDataset object also contains functions that apply operations to the entire dataset, or that write output for a range of samples/runs.

Building a FacetsDataset object involves three steps:
* Create a FacetsDataset from a built FacetsMeta: `test_dataset = FacetsDataset(prepared_metadata)`
* Apply desired filter/selection functions.
* Execute `buildFacetsDataset()`.

Filters that need clinical data (cancer type, detail, clinical purity, onko code, TMB, MSI) raise `FacetsConfigError` when the FacetsMeta was built without a clinical sample file. Clinical values are looked up by tumor sample id whatever the run key shape (single-run, all-fits, 2n classes).

#### Selection / Filtration Functions
* setCancerTypeFilter(string[] selectedCancerTypes)
  * This function applies a filter to the FacetsDataset build so that only indicated cancer types will be selected. Input is expected to a list of strings, even if the length of the selection is 1. Note that cancer types are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setCancerTypeFilter(["Example Cancer Type 1", "Example Cancer Type 2"])`.
* setCancerTypeDetailFilter(string[] selectedCancerDetailTypes)
  * This function applies a filter to the FacetsDataset build so that only indicated detail cancer types will be selected. Input is expected to a list of strings, even if the length of the selection is 1. Note that cancer types are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setCancerTypeDetailFilter(["Example Cancer Type 1", "Example Cancer Type 2"])`.
* setPurityFilter(minPurity, maxPurity=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples with purity values in a given range are selected.  minPurity is required, however maxPurity, unless otherwise indicated, defaults to 1. Runs whose FACETS purity is exactly 0.3 (FACETS' value when it could not estimate purity) or NA are excluded.
* setClinicalPurityFilter(minPurity, maxPurity=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples with clinical purity values in a given range are selected. Clinical purity is the purity value extracted from the data_clinical_sample.txt file, whereas purity, as used in setPurityFilter(), references FACETS calculated purity. minPurity is required, however maxPurity, unless otherwise indicated, defaults to 1. Samples with NA clinical purity are excluded.
* setOnkoCodeFilter(selectedOnkoCodes)
  *  This function applies a filter to the FacetsDataset build so that samples with specified onko codes are included. Input is expected to a list of strings, even if the length of the selection is 1. Note that onko codes are case sensitive and should correspond to values in the data_clinical_sample.txt file. I.E. `setOnkoCodeFilter(["LUAD","LUSC"])`.
* setPloidyFilter(minPloidy, maxPloidy)
  *  This function applies a filter to the FacetsDataset build so that only samples within a given ploidy range are selected. minPloidy and maxPloidy are both required.
* setDipLogRFilter(minDipLogR, maxDipLogR)
  *  This function applies a filter to the FacetsDataset build so that only samples within a given dipLogR range are selected. minDipLogR and maxDipLogR are both required.
* setCvalFilter(minC, maxC)
  * This function applies a filter to the FacetsDataset build so that only samples within a given cValue range are selected. minC and maxC are both required.
* setWGDFilter(keepWGD)
  * This function applies a filter to the FacetsDataset build indicating whether to select only WGD samples, or only non-WGD samples.  If this is set to True, WGD samples will be selected.  If set to False, only non-WGD samples will be selected.
* setFGAFilter(minFGA, maxFGA=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples within a given FGA range are selected.  minFGA is required, however maxFGA, unless otherwise indicated, defaults to 1.
* setFracLohFilter(minLoh, maxLoh=1.0)
  * This function applies a filter to the FacetsDataset build so that only samples within a given frac_loh range are selected.  minLoh is required, however maxLoh, unless otherwise indicated, defaults to 1.
* setFacetsQCFilter(keepQCPass)
  * This function applies a filter to the FacetsDataset build indicating whether to select only samples that pass FACETS QC, or only those who fail.  If this is set to True, only FACETS QC PASS samples will be selected. If set to False, only FACETS QC FAIL samples will be selected. The flag is read from the sample's `facets_qc.txt` row for the run's fit (nested ultra names such as `default/ultra` are matched as written).
* setTMBFilter(minTMB, maxTMB)
  * This function applies a filter to the FacetsDataset build so that only samples within a given TMB range are selected. minTMB and maxTMB are both required.
* setMSIFilter(minMSI, maxMSI)
  * This function applies a filter to the FacetsDataset build so that only samples within a given MSI range are selected. minMSI and maxMSI are both required.


#### Dataset Functions
* buildFacetsDataset()
  * This function will populate the FacetsDataset with samples and/or runs that meet the specified selection criteria indicated by using the above filtration functions, using the directory mappings in the [FacetsMeta](facetsmeta.md) object it was created with.  Two structures hold the result: `sampleList`, a dictionary `{pair_id -> FacetsSample}` keyed by the `<tumor>_<normal>` pair id, and `runList`, a list of `FacetsRun`.  A [FacetsSample](facetssample.md) groups the runs of one pair (alt fits in all-fits mode, or the clinical and research runs of a 2n pair with `fit_class="both"`), whereas a [FacetsRun](facetsrun.md) represents a single fit. Both reference the same objects. After the build, `numSamples` (runs), `num_filtered`, `num_file_failed` and `num_build_failed` summarize what happened; a run whose files are malformed is skipped and recorded in `FacetsMeta.failed_samples`.
* runAlterationAnalysis()
  * This function will run alteration analysis on a constructed FacetsDataset object.  This process will run several functions in each [FacetsRun](facetsrun.md).  Specifically, for each sample, arms will be defined for all [FacetsSegments](facetssegment.md), then LoH and Gains will be calculated at the arm level, CF levels will be calculated, and run level alteration calls will be made.  (i.e. "Hypoploidy")
* writeAlterationData(outfile_path)
  * This function will write alteration data for all samples in the FacetsDataset to the specified outfile_path.  This includes arm level gain/loss calls, CF levels, and run-level calls.
* writeReport(outfile_name)
  * Writes one row per run: ID, fitDir, fitType (`Best`, `Acceptable`, `AutoQC_Best`, `AutoQC_Pass`, `Default`, `Unreviewed`), clinical and FACETS metrics.
* writePurityCFs(outfile_path, use_base_cf=False, use_long_id=True)
  * This function will calculate and write calculated purity data for each sample.  Purity is defined as the greatest CF value that is not NA or 1. If use_base_cf is true, cf will be used rather than cf.em.  If use_long_id is True, the pair id including the paired normal is written instead of the run key.
* printFacetsSampleById(id)
  * This function will print a specific sample in this facets dataset. id is the pair id (`<tumor>_<normal>`).
* writeDatasetSummary(outdir=".")
  * This function will output a text file containing information that summarizes a facets dataset. This includes metrics like number of best/acceptable fits, number of failed samples that couldn't build, number of cancer types, etc.  The output's name contains the date/time; it is written into `outdir`.
* createHistogram(variable, inbins=10, outdir=".")
  * This function creates a histogram png from the requested variable and returns its path.  It defaults to 10 automated bins for numerical variables but that can be changed. For strings like cancerType or boolean values like WGD it uses however many bins are needed.
  * Available variables are as follows: purity, clinicalPurity, ploidy, dipLogR, cval, fga, frac_loh, tmb, msi, cancerType, onkoCode (oncoCode accepted), wgd, review_status, fit_class, fit_name.
* copyDatasetToFolder(output_folder)
  * Copies every sample directory of the dataset to `<output_folder>/all/<P-00xxx>/<pair>/` without overwriting (the whole pair directory for 2n).
* backupMafs(ref_meta)
  * Date-stamped `.bak` copies of every run's ccf.maf.


#### Examples

```python
  # Build a FacetsDataset that contains only Breast Cancer type samples with purity between 0.2 and 0.8.
   prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
   prepared_metadata.setSingleRunPerSample(True, True)
   prepared_metadata.buildFacetsMeta()
   test_dataset = FacetsDataset(prepared_metadata)
   test_dataset.setPurityFilter(0.2, 0.8)
   test_dataset.setCancerTypeFilter(["Breast Cancer"])
   test_dataset.buildFacetsDataset()
```

```python
  # Build a FacetsDataset over ALL fits of a few onkocode specific samples with purity between 0.3 and 1,
  # and run alteration analysis. Output the result to file.
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/path/to/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "hisens")
  prepared_metadata.buildFacetsMeta()
  test_dataset = FacetsDataset(prepared_metadata)
  test_dataset.setPurityFilter(0.3)
  test_dataset.setOnkoCodeFilter(["SCLC","CSCLC"])
  test_dataset.buildFacetsDataset()
  test_dataset.printFacetsSampleById("P-0000004-T01-IM3_P-0000004-N01-IM3")
  test_dataset.runAlterationAnalysis()
  test_dataset.writeAlterationData("test_alterations_data.txt")
```

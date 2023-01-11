## Integrity Check

### 
This script checks the integrity of an inputted facetsMeta object or facetsDataset object.  

integrity_check(inFacetsMeta = False, inFacetsDataset = False) 

### Inputs: 
* inFacetsMeta (FacetsMeta object):  Inputted facets meta object that has been constructed previously.  If no inFacetsDataset is input then the integrity check will check the facets directory indicated in FacetsMeta 
    
* inFacetsDataset (facetsDataset object) : Inputted facetsDataset object that has been constructed previously. The value inFacetsDataset.ref_facetsMeta will replace te input FacetsMeta object if there was one. 

One of these must be entered as an arguement.  If entering a facetsDataset then FacetsMeta information will be gathered from the facetsDataset as needed.  If entering a FacetsMeta object then the integrity of the facets repo associated with that facetsMeta object will be the subject of the integrity check.  

### Outputs
FACETSintegrity_{datetime}.tsv : A TSV containing the data collected from integrity.py

Cols are as follows : 
* ID : Name of the sample folder
* emptyfiles:  Files that have a size of 0 bytes but do exist
* missing_filesin_samplefolder: Files that do not exist in the sample folder and should be there
* MissingfilesinFitfolder:  Files that do not exist in the fit folder but should be there
* missingcols:  Files that are missing columns from the fit folder
* tumorbam_size_pass: The tumor bam file is within two standard deviations of the mean file size
* normalbam_size_pass: The tumor bam file is within two standard deviations of the mean file size
* manifest_file_labelerror: Manifest files that have multiple fits of the same type
* manifest_folder_missing:  A fit in the manifest does not have the fit folder
* Consent:  12245-a consent.  1 if consented
* FitCount: Number of fits
* In_Clinical_file:  If constructed from facetsMeta this will be 1 if the sample is in the clincal file. 

Also in terminal are printed several metrics such as: 
* Total number of samples
* Number of fits
* Number of Best, acceptable, default, and no fit
* Number of samples with multiple normals but same tumor 


### Usage Example
```python
def test_Integritycheck(useSingleRun, allowDefaults,clinsample, facets_dir ):
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
    #facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity")
    prepared_metadata.setSingleRunPerSample(useSingleRun,allowDefaults)
    prepared_metadata.buildFacetsMeta()

    #check Integrity
    integrity_check(inFacetsMeta = prepared_metadata)

    ## To test the dataset function of integrity check
    #test_dataset = FacetsDataset(prepared_metadata)
    #test_dataset.setPurityFilter(.5,1)
    #test_dataset.buildFacetsDataset()

    ## Check Integrity
    # integrity_check(inFacetsDataset = test_dataset)


clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
facets_dir            = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/run_11-14-22/facets/fix_alt/all/"
#facets_dir           = "/work/ccs/shared/resources/impact/facets/all/"
test_Integritycheck(True, True, clinical_sample_file, facets_dir)
```

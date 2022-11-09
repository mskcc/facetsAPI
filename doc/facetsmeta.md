## FacetsMeta

Most tasks using facetsAPI will begin with creating a FacetsMeta object that represents the FACETS repository you are working with. 
Primarily, this class encapsulates clinical data and maps the target FACETS repository directory structure to identify and store paths to relevant 
files and folders in the FACETS repository.

### Initializing FacetsMeta

The initialization function of a FacetsMeta class requires 4 parameters:
* Data clinical sample file path.
* FACETS repository base path.
* Use default fit when best/acceptable fits not available (Boolean)
* Use "purity" or "hisens" data.
* A path to store a persistent data file for this object. (Optional)

The **data clinical sample file** is updated daily by the DMP-2022 knowledgebase and can be obtained through their github page.  
Note that not all samples present in the data_clinical_sample.txt file will necessarily be present in a FACETS file repository depending on
factors such as change in consent, most recent FACETS update run date, or custom repositories.  

The **FACETS repository base path** represents the base directory of a FACETS repository file system. 
For example, the main Impact FACETS repository is located on Juno at `/work/ccs/shared/resources/impact/facets/all/`.

The **Use default fit when best/acceptable fits not available** is a Boolean parameter (True/False).  
If this is set to `False`, the facets_manifest.txt file will be checked for each sample and only fits with "Best Fit" or "Accepted Fit" will be included
in subsequent data processing steps. Default fits will not be included unless they are also Best or Accepted.  If this is `True`, Best or Accepted fits
will be chosen as a priority when building data sets, but in the absense of an ideal fit, the default fit will be included in this FacetsMeta object.

The **Use "purity" or "hisens" data** parameter will determine which type of data is examined when looking at run level data.  
In the primary FACETS repository, hisens runs (cVal=50) and purity runs (cVal=100) both exist with *_purity or *_hisens file name suffixes.  
This parameter allows you to decide which of these run types to select for this FacetsMeta object.  Valid selections are "purity" or "hisens".

The **path to store a persistent data file** parameter is optional.  
If no value is provided for this parameter, the data_clinical_sample.txt file will be read and the indicated FACETS repository will be scanned and 
structured into a FacetsMeta object.  This process can take a few minutes, and it is sometimes not necessary to repeat the process for each run.  
By providing a path to this parameter, on the initial run, a binary version of the fully structured FacetsMeta object will be stored to the indicated path.
In subsequent runs, if the indicated path already points to a file, instead of parsing the clinical data and scanning the facets repository, the 
previously processed data can simply be loaded into the FacetsMeta object.  


#### Examples

```python
  # This will build a FacetsMeta object that includes default fits, selects "purity" type data, 
  # and persists the FacetsMeta object into a file called all_meta.dat
  
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, True, "purity", "all_meta.dat")
```


```python
  # This will build a FacetsMeta object that includes only curated fits, 
  # selects "hisens" type data, and does not persist data.
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, False, "hisens")
```

### Functions

* printLogo() 
  * A simple method for printing a FACETS API logo. 
  **This function is run automatically during FacetsMeta construction**.
* parseClinicalSample()
  * This function will accept a clinical sample file and scan each sample's corresponding facets directory.
  the facets_review.manifest file will be read to identify
  the appropriate fit to use in analysis.  A dictionary mapping
  samples to their appropriate directories is built in master_file_dict.
  This function will also build dictionaries mapping sample id to
  cancer type and to patient id. 
  **This function is run automatically during FacetsMeta construction**.




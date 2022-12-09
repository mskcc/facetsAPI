## FacetsMeta

Most tasks using facetsAPI will begin with creating a FacetsMeta object that represents the FACETS repository you are working with. 
Primarily, this class encapsulates clinical data and maps the target FACETS repository directory structure to identify and store paths to relevant 
files and folders in the FACETS repository.

### Initializing FacetsMeta

The initialization function of a FacetsMeta class requires 4 parameters:
* Data clinical sample file path. (Optional)
* FACETS repository base path.
* Use "purity" or "hisens" data.
* A path to store a persistent data file for this object. (Optional)

The **data clinical sample file** is updated daily by the DMP-2022 knowledgebase and can be obtained through their github page.  
Note that not all samples present in the data_clinical_sample.txt file will necessarily be present in a FACETS file repository depending on
factors such as change in consent, most recent FACETS update run date, or custom repositories.  **NOTE: If this parameter is an empty string, then the provided FACETS repository base path will be scanned to produce the set of samples present in the directory.  This will mean that all samples/fits in the folder will be included in the FacetsMeta object, but clinical metadata from the data_clinical_sample file will not be available, and some downstream functions will be restricted.**

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

### Functions

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
* setSingleRunPerSample(bool useSingleSample, bool allowDefaults)
  * This function will set this FacetsMeta object to select a single run per sample based on identified best or acceptable fits, as defined in the facets manifest for any given sample.  If the first parameter, useSingleSample, is set to True, then for each sample, a single run will be selected if a best or acceptable fit is selected. If useSingleSample is set to True, then it is also possible to allow default fits to be selected in cased where no best or acceptable fit is indicated. Default for both of these parameters is set to False.
* setVerbose(bool activateVerbose)
  * This function accepts True or False as a parameter, and will activate verbose mode if set to True.  This will make warning messages visible during runtime and show more details of data processing for various processes.  By default, verbose mode is set to False.   


### Additional Data

* OnkoTree Code generalized lists are included and can be referenced by the following variables.  For example `FacetsMeta.nscLung_cancer`.

```python
    #OncoTree Codes that correspond to a specific cancer type.  
    breast_carcinoma    = ["ILC","IDC","BRCA","BRCNOS","BRCANOS","MDLC","MBC","CSNOS"]
    nscLung_cancer      = ["LUAD","LUSC","LUNE","NSCLC","LUAS","NSCLCPD","ALUCA","SARCL"]
    endometrial_cancer  = ["UMEC","UCCC","UEC","UCEC","UCS","UDDC","UUC","USC","OUTT"]
    colorectal_cancer   = ["READ","COAD","MACR","COADREAD","CMC","CM"]
    melanoma            = ["SKCM","ACRM","MUP","ARMM","HNMUCM","UM","VMM","SKCN"]
    scLung_cancer       = ["SCLC","CSCLC"]
    ovarian_cancer      = ["HGSOC","CCOV","LGSOC","OCS","EOV"]
    prostate_cancer     = ["PRAD","PRSCC","PRSC"]
    unknown_primary     = ["NETNOS","ADNOS","CUP","PDC","NECNOS","CUPNOS","SCUP"]
    bladder_cancer      = ["BLCA", "UTUC", "UCU", "BLSC", "USCC"]
    pancreatic_cancer   = ["PAAD","PANET","PAASC","UCP","SPN","PB"]
    renalCell_carcinoma = ["PRCC","TRCC","CCRCC","URCC","RCC","CHRCC","MT","SRCC","ROCY","MTSCC"]
    glioma              = ["ASTR","ODG","AODG","GBM","HGGNOS","AASTR","GB","DIFG"]
    headNeck_carcinoma  = ["HNSC","OPHSC","HPHSC","OCSC","HNSCUP","LXSC","HNNE","SNSC","ODGC"]
    germCell_tumor      = ["MGCT","SEM","VMT","BMT","OYST","GCTSTM","NSGCT","EMBCA","OGCT","OMGCT","TT","TYST","VDYS","ODYS","OIMT","VYST","VMGCT","BMGCT","BIMT","VIMT","BYST","OMT","GCT"]

```


* chr_arms
  * Chromosome arm position map.  This structure is a dictionary mapping chromosome to relevant arm level positions. chrID -> [p_start, p_end, q_start, q_end].  For example, chromosome 1 data is  `{1: [1,125000000,125000001,249250621]}`.  Data for chromosome 1, for example, can be accessed using `FacetsMeta.chr_arms.get(1)`. 


#### Examples

```python
  # This will build a FacetsMeta object, selecting a best/acceptable run for each sample 
  # and includes default fits when best/acceptable is not indicated. Selects "purity" type data, 
  # and persists the FacetsMeta object into a file called all_meta.dat
  
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity", "all_meta.dat")
  prepared_metadata.setSingleRunPerSample(True, True)
  prepared_metadata.buildFacetsMeta()
  
```

```python
  # This will build a FacetsMeta object, selecting a best/acceptable run for each sample 
  # and includes default fits when best/acceptable is not indicated. With the clinical
  # sample file not provided, the run will be based on what is present in the facets_dir.
  # Note that some clinical metadata will be unavailable in this scenario for downstream analysis.
  
  clinical_sample_file  = ""
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity", "all_meta.dat")
  prepared_metadata.setSingleRunPerSample(True, True)
  prepared_metadata.buildFacetsMeta()
  
```

```python
  # This will build a FacetsMeta object that includes all samples and runs. 
  # selects "hisens" type data, and does not persist data.
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "hisens")
  prepared_metadata.buildFacetsMeta()
```

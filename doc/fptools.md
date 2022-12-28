## FPTools

The FPTools class functions as a wrapper to execute certain functions from Facets Preview and execute some types of file operations on a facets directory. 

### Initializing a FPTools Object

A FPTools object can be instantiated by providing the following variables:

* fp_config - A path to a facets preview config file. I.E. "/path/to/facets_preview_conf.json"
* cbio_maf - A path to a vep88 maf file. I.E. "/path/to/data_mutations_extended.oncokb.vep.maf"
* cbio_nonsigned_maf - A path to a vep88 nonsignedout maf file. I.E. "/path/to/data_mutations_nonsignedout.vep.maf"

An example of an instantiation command for a FPTools object is:
`fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf)`

### Functions

* loadModule
  * This function will run a system command to load a specific module on Juno.  For example to load R-3.6.3, `fp_tools.loadModule("R/R-3.6.3")`
* runGenerateGenomicAnnotations(short_id, long_id, path)
  * This function will run Facets Preview's generate genomic annotations function on a specified sample.  Short ID is a sample's DMP id with only tumor (P-012345-T01), long ID includes the normal (P-012345-T01_P-012345-N01), path is the base path to the facets sample directory.  An example of this function's use can be seen in [this script](../scripts/correct_missing_anno.py), which will attempt to run this function on any samples missing genomic annotations.  This function will submit LSF jobs to the juno cluster for the requested sample run.
* backupManifests(ref_meta)
  * This function will make backup files of existing manifests for every sample in the provided FacetsMeta object with a date stamped file name.
* backupFacetsQC(ref_meta)
  * This function will make backup files of existing facets qc files for every sample in the provided FacetsMeta object with a date stamped file name.


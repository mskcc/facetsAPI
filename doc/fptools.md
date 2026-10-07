## FPTools

The FPTools class functions as a wrapper to execute certain functions from Facets Preview and execute some types of file operations on a facets directory.

### Initializing a FPTools Object

A FPTools object can be instantiated by providing the following variables:

* fp_config - A path to a facets preview config file. I.E. "/path/to/facets_preview_conf.json". May also come from `FACETSAPI_FP_CONFIG` or `FacetsConfig.fp_config_file`.
* cbio_maf - A path to a vep88 maf file. I.E. "/path/to/data_mutations_extended.oncokb.vep.maf"
* cbio_nonsigned_maf - A path to a vep88 nonsignedout maf file. I.E. "/path/to/data_mutations_nonsignedout.vep.maf"
* config (keyword, optional) - a [FacetsConfig](configuration.md) controlling how `Rscript` is reached and where the R wrapper and the facetsPreview R library live.

```python
fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf)
fp_tools = FPTools(fp_config, cbio_maf, cbio_nonsigned_maf,
                   config=FacetsConfig(r_prefix="module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3",
                                       r_lib_loc="/path/to/R_libs"))
```

### Functions

* runGenerateGenomicAnnotations(short_id, long_id, path)
  * This function will run Facets Preview's generate genomic annotations function on a specified sample. Short ID is a sample's DMP id with only tumor (P-012345-T01), long ID includes the normal (P-012345-T01_P-012345-N01), path is the base path to the facets sample directory (for a 2n pair, the class directory). An example of this function's use can be seen in [this script](../scripts/correct_missing_anno.py). The R wrapper `R/fp_generate_genomic_anno_wrapper.R` runs inline through `bash -lc`; a nonzero exit raises `ExternalToolError`.
* backupManifests(ref_meta)
  * This function will make backup files of existing manifests for every sample in the provided FacetsMeta object with a date stamped file name. Returns the number of files backed up.
* backupFacetsQC(ref_meta)
  * This function will make backup files of existing facets qc files for every sample in the provided FacetsMeta object with a date stamped file name. Returns the number of files backed up.

`loadModule()` was removed in 0.6.0 (see [ExtTools](exttools.md)).

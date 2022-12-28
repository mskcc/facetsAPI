## ExtTools

The ExtTools class functions as a wrapper to execute certain functions from misc. external sources. 

### Initializing a ExtTools Object

A ExtTools object can be instantiated by providing a built [FacetsMeta](facetsmeta.md) object.

An example of an instantiation command for a ExtTools object is:
`ext_tools = ExtTools(prepared_metadata)`

### Functions

* loadModule
  * This function will run a system command to load a specific module on Juno.  For example to load R-3.6.3, `fp_tools.loadModule("R/R-3.6.3")`
* makeMergedFile(fileType, outfile_path)
  * This function will make a merged file for any file type in the MetaDictMap with a shared header. fileType should be a [MetaDictMap](metadictmap).FILETYPE. This will write the merged file to the designated outfile_path and return a path to the written file if successful.
* runAscets(output_dir="", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True)
  * This function will run ascets on a sample by submitting the LSF job to the Juno cluster.  Ascets files required for this function are stored in the [tools/ascets directory](../tools/ascets/) of this repository.  Default values for required parameters are provided in line with guidance from the [Ascets documentation](https://github.com/beroukhim-lab/ascets). use_adjusted is a toggle that configures whether to use adjusted or unadjusted .seg file data for processing.  



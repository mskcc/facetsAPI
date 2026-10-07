## ExtTools

The ExtTools class functions as a wrapper to execute certain functions from misc. external sources.

### Initializing a ExtTools Object

A ExtTools object can be instantiated by providing a built [FacetsMeta](facetsmeta.md) object and, optionally,
a [FacetsConfig](configuration.md) describing how to reach `Rscript`:

```python
ext_tools = ExtTools(prepared_metadata)                                   # FacetsConfig.from_env()
ext_tools = ExtTools(prepared_metadata, config=FacetsConfig(r_prefix="module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3"))
```

R calls run **inline** (no job scheduler). Wrap your script in `sbatch`, a Nextflow process or a container if you need one.

### Functions

* makeMergedFile(fileType, outfile_path)
  * This function will make a merged file for any file type in the MetaDictMap with a shared header. fileType should be a [MetaDictMap](metadictmap.md).FILETYPE. This will write the merged file to the designated outfile_path and return a path to the written file if successful. Files shared by several runs of the same sample (alt fits, 2n classes) are merged once.
* runAscets(output_dir="", ref_genome_coords="hg19", min_arm_breadth=0.5, keep_noise=False, arm_alt_frac_thresh=0.7, use_adjusted=True)
  * This function will run ascets over the merged seg file of this object's FacetsMeta. Ascets files required for this function are stored in the [tools/ascets directory](../tools/ascets/) of this repository (override with `FacetsConfig.ascets_dir` / `FACETSAPI_ASCETS_DIR`). Default values for required parameters are provided in line with guidance from the [Ascets documentation](https://github.com/beroukhim-lab/ascets). use_adjusted is a toggle that configures whether to use adjusted or unadjusted .seg file data for processing. `runLog.out` / `runLog.err` are written to output_dir; a nonzero R exit raises `ExternalToolError`. Returns the output prefix.
* runFacets()
  * Not implemented (raises `NotImplementedError`). FACETS runs are produced by the IMPACT FACETS pipeline, which this API only reads.

`loadModule()` was removed in 0.6.0: `os.system("module load ...")` ran in a throwaway shell and never affected the later `Rscript` call. Use `FacetsConfig.r_prefix` instead.

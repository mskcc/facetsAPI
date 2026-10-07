## Configuration (FacetsConfig)

facetsAPI itself is pure Python. Only `ExtTools.runAscets()` and `FPTools.runGenerateGenomicAnnotations()` call R,
and where R lives differs per system. `FacetsConfig` describes how to reach it.

| setting | env var | default | meaning |
|---|---|---|---|
| `r_prefix` | `FACETSAPI_R_PREFIX` | `""` | shell snippet run **in the same shell** before `Rscript`, e.g. `module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3` |
| `rscript` | `FACETSAPI_RSCRIPT` | `Rscript` | the Rscript executable |
| `fp_wrapper_script` | `FACETSAPI_FP_WRAPPER` | `<repo>/R/fp_generate_genomic_anno_wrapper.R` | facets-preview annotation wrapper |
| `fp_config_file` | `FACETSAPI_FP_CONFIG` | none | facets-preview config json (required by FPTools) |
| `r_lib_loc` | `FACETSAPI_R_LIB_LOC` | none | R library directory holding facetsPreview; passed to the wrapper |
| `ascets_dir` | `FACETSAPI_ASCETS_DIR` | `<repo>/tools/ascets` | `run_ascets.R` and the genomic arm coordinate files |

Precedence: explicit constructor argument > environment variable > default.

```python
from facetsAPI import FacetsConfig, ExtTools

cfg = FacetsConfig.from_env()                                   # env vars + defaults
cfg = FacetsConfig.from_env(rscript="/opt/R/4.4.1/bin/Rscript")  # explicit override wins
print(cfg.compose_r_command(["/x/run_ascets.R", "-i", "merged.seg"]))
# module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3 && Rscript /x/run_ascets.R -i merged.seg
```

Every R call is composed into one string and executed with `bash -lc`, so the prefix and `Rscript` share an
environment. On IRIS, Lmod's `module` is an exported shell function and is available there. A nonzero exit raises
`ExternalToolError`, which carries the command, the return code and the tail of stderr.

Typical IRIS shell setup:

```bash
export FACETSAPI_R_PREFIX="module use /admin/software/spack/spack_modulefiles/Core/ && module load r/3.6.3"
export FACETSAPI_FP_CONFIG=/path/to/config_facets_preview.json
export FACETSAPI_R_LIB_LOC=/path/to/R_libs
```

The defaults for the wrapper script and the ascets directory point into this repository checkout, so an editable
install (`pip install -e .`) or a `PYTHONPATH` checkout works out of the box; for other install layouts set
`FACETSAPI_FP_WRAPPER` and `FACETSAPI_ASCETS_DIR`.

### Exceptions

All errors raised by the library derive from `FacetsAPIError`:

| class | raised for |
|---|---|
| `FacetsConfigError` | bad arguments (invalid `fit_class`, `hisens_vs_purity`, filter arguments, missing config file) |
| `FacetsRepositoryError` | repository path or clinical sample file missing/unreadable/lacking columns |
| `ManifestError` (a `FacetsRepositoryError`) | a `facets_review.manifest` is unreadable or lacks required columns |
| `ParseError` | a FACETS output file we chose to load has malformed content |
| `ExternalToolError` | an R call exited nonzero |

Per-sample problems (no manifest, a fit missing files, a malformed fit) are **not** exceptions: they are counted in
`FacetsMeta.build_counts`, listed in `FacetsMeta.failed_samples`, printed when `setVerbose(True)`, and the sample is skipped.

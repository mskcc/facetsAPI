## MetaDictMap

This class is a helper class to keep track of the types and order of files in a [FacetsMeta](facetsmeta.md) master_file_dict type object.  
As the master_file_dict keeps track of a large number of files for each sample, it is complicated to reference these values by index.  
When calling master_file_dict, instead of numbered indices, use MetaDictMap.FILE_NAME.

```
OUT_FILE = 0             # <pair>_purity.out (or _hisens)
CNCF_FILE = 1
QC_FILE = 2              # <pair>.qc.txt (fit level)
FACETS_QC_FILE = 3       # facets_qc.txt (sample / 2n class level)
SELECTED_FIT_DIR = 4     # the fit directory, resolved from the discovered sample directory (trailing slash)
GENE_FILE = 5
ADJUSTED_SEG_FILE = 6
CCF_MAF_FILE = 7
NONSIGNED_MAF_FILE = 8
UNADJUSTED_SEG_FILE = 9
MANIFEST_FILE = 10
FIT_STATUS = 11          # "review_status#fit_name#reviewed_by"
SAMPLE_BASE_DIR = 12     # the directory holding the manifest: the sample dir (standard) or the class dir (2n)
RDATA_FILE = 13          # "" when the fit has no .Rdata
FIT_CLASS = 14           # "standard" | "clinical" | "research"           (added in 0.6.0)
FIT_TIER = 15            # "reviewed" | "autoqc" | "default" | "unreviewed" (added in 0.6.0)
```

`FIT_TIER` is the tier the fit was selected under in single-run mode (see [fit selection](facetsmeta.md#fit-selection-policy)),
or the row's natural tier in all-fits mode.

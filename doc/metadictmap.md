## MetaDictMap

This class is a helper class to keep track of the types and order of files in a [FacetsMeta](facetsmeta.md) master_file_dict type object.  
As the master_file_dict keeps track of a large number of files for each sample, it is complicated to reference these values by index.  
When calling master_file_dict, instead of numbered indeces, use MetaDictMap.FILE_NAME.

```
OUT_FILE = 0
CNCF_FILE = 1
QC_FILE = 2
FACETS_QC_FILE = 3
SELECTED_FIT_DIR = 4
GENE_FILE = 5
ADJUSTED_SEG_FILE = 6
CCF_MAF_FILE = 7
NONSIGNED_MAF_FILE = 8
UNADJUSTED_SEG_FILE = 9
MANIFEST_FILE = 10
FIT_STATUS = 11
SAMPLE_BASE_DIR = 12
```

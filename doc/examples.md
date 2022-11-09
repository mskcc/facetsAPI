## Examples

#### Selecting a FACETS data subset.

The following code will build a FacetsDataset with a set of filters applied so that only samples passing FACETS QC, 
who have been specified as "best fits", and with a minimum purity cutoff of 0.2.  After selecting and building this dataset the user could
manipulate or output it as they require.

```
from facetsAPI import *

if __name__ == '__main__':
    clinical_sample_file  = "/work/ccs/shared/resources/impact/cbio_mutations/adam_cron/bsub_run/data_clinical_sample.oncokb.txt"
    facets_dir            = "/work/ccs/shared/resources/impact/facets/all/"

    #Initialize FacetsMeta. This will build all relevant metadata we need going forward.
    prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, True, "purity", "all_meta.dat")

    hypo_dataset = FacetsDataset()
    filter_on_qc_fail       = True     # If this is True, facets_qc failed samples will be excluded from analysis.
    purity_filter_thresh    = 0.2      # This is the level of purity at which we disqualify samples from consideration.
    hypo_dataset.setPurityFilter(purity_filter_thresh)
    hypo_dataset.setFacetsQCFilter(filter_on_qc_fail)
    hypo_dataset.buildFacetsDataset(prepared_metadata)

```

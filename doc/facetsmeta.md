## FacetsMeta

Most tasks using facetsAPI will begin with creating a FacetsMeta object that represents the FACETS repository you are working with.
Primarily, this class encapsulates clinical data and maps the target FACETS repository directory structure to identify and store paths to relevant
files and folders in the FACETS repository.

### Initializing FacetsMeta

```python
FacetsMeta(clinical_sample_file="", facets_repo_path="", hisens_vs_purity="purity", persist_data="no",
           fit_class="research", fit_policy=None)
```

The **data clinical sample file** is updated daily by the DMP-2022 knowledgebase and can be obtained through their github page.
Note that not all samples present in the data_clinical_sample.txt file will necessarily be present in a FACETS file repository depending on
factors such as change in consent, most recent FACETS update run date, or custom repositories.  **NOTE: If this parameter is an empty string, then the provided FACETS repository base path will be scanned to produce the set of samples present in the directory.  This will mean that all samples/fits in the folder will be included in the FacetsMeta object, but clinical metadata from the data_clinical_sample file will not be available, and some downstream functions will be restricted.**

The **FACETS repository base path** represents the base directory of a FACETS repository file system: the directory that holds the
`P-00xxx` patient folders. The repository can be mounted anywhere; nothing in facetsAPI depends on the cluster path recorded inside
the manifests.

The **hisens_vs_purity** parameter will determine which type of data is examined when looking at run level data.
In the primary FACETS repository, hisens runs (cVal=50) and purity runs (cVal=100) both exist with *_purity or *_hisens file name suffixes.
Valid selections are "purity" or "hisens".

The **persist_data** parameter is optional.
If no value is provided for this parameter, the data_clinical_sample.txt file will be read and the indicated FACETS repository will be scanned and
structured into a FacetsMeta object.  This process can take a few minutes, and it is sometimes not necessary to repeat the process for each run.
By providing a path to this parameter, on the initial run, a binary version of the fully structured FacetsMeta object will be stored to the indicated path.
In subsequent runs, if the indicated path already points to a file, the previously processed data is loaded instead. Files written by facetsAPI 0.5
load with a warning; the fit-selection settings then come from the constructor and setters rather than the file.

The **fit_class** parameter only matters for **2n repositories** (see below): `"clinical"`, `"research"` (default) or `"both"`.

The **fit_policy** parameter is a [FitPolicy](#fit-selection-policy) controlling which fit is chosen per sample in single-run mode.

### Repository layouts: standard and 2n

A **standard** sample directory holds `facets_review.manifest`, `facets_qc.txt` and one directory per fit:

```
<repo>/P-00000/P-0000001-T01-IM6_P-0000001-N01-IM6/facets_review.manifest
                                                  /facets_qc.txt
                                                  /default/  alt_dipLogR_-0.03/  refit_.../
```

A **2n** pair (the facets2n / impact_2n repository) has no root manifest. Its `clinical/` and `research/` subdirectories are each a
complete standard-shaped sample directory, and research fits may carry an `ultra` sub-fit (`default/ultra`, hisens-only):

```
<repo>/P-00000/<pair>/clinical/{facets_review.manifest, facets_qc.txt, default/, alt_dipLogR_-0.26/}
                     /research/{facets_review.manifest, facets_qc.txt, default/, default/ultra/, alt_dipLogR_-0.26/, ...}
```

The layout is detected per sample, so a FacetsMeta works on either repository without configuration. For 2n pairs `fit_class`
decides which class is loaded; `"both"` loads each class as its own run. Run keys and ids:

| mode | standard sample | 2n pair, one class | 2n pair, `fit_class="both"` |
|---|---|---|---|
| single run | `P-0000001-T01-IM6` | `P-0000007-T01-IM6` | `P-0000007-T01-IM6#clinical`, `P-0000007-T01-IM6#research` |
| all fits | `<pair>#default` | `<pair>#research/default` | `<pair>#clinical/default`, `<pair>#research/default` ... |

`FacetsMeta.shortId(key)` reduces any key to the tumor sample id; `pairIdForKey(key)` gives the `<tumor>_<normal>` pair id that
[FacetsDataset](facetsdataset.md) groups runs under. Each run records its class at `MetaDictMap.FIT_CLASS` and on `FacetsRun.fit_class`.
Ultra fits are excluded by default (`FitPolicy(include_ultra_fits=True)` admits them; note that on disk they carry only `_hisens` files and no
qc / gene-level files, so they load only when those exist).

### Fit selection policy

Single-run mode picks one fit per sample (per class for 2n) from the manifest, in tiers:

| tier | manifest `review_status` | meaning |
|---|---|---|
| `reviewed` | `reviewed_best_fit` then `reviewed_acceptable_fit` | a human reviewer chose the fit |
| `autoqc` | `auto_qc_best_fit` then `auto_qc_pass` | the autoQC ruleset chose the fit |
| `default` | any status, `fit_name == "default"` | facets' own first fit |

The default policy is `FitPolicy(("reviewed", "autoqc", "default"))`. Selection **falls through**: the first candidate whose files
exist for the requested purity/hisens mode wins, so a sample whose only autoQC row is an ultra fit still gets its default fit in purity mode.
Rows with `fit_name == "Not selected"` carry no fit and are ignored.

```python
from facetsAPI import FacetsMeta, FitPolicy

m = FacetsMeta(clin, repo, "purity")
m.setSingleRunPerSample(True, True)                # reviewed > autoqc > default
m.setSingleRunPerSample(True, False)               # reviewed > autoqc (no defaults)
m.setFitPolicy(FitPolicy.reviewed_only())          # only human-reviewed fits
m.setFitPolicy(FitPolicy.autoqc_only())            # only autoQC fits
m.setFitPolicy(FitPolicy.default_only())           # only the default fit
m.setFitPolicy(FitPolicy(("autoqc", "reviewed", "default"), tie_break="newest_reviewed"))
m.setFitPolicy(FitPolicy(include_ultra_fits=True))
```

A policy passed to the constructor or `setFitPolicy()` survives a later `setSingleRunPerSample(True, allowDefaults)`; only its default tier is toggled.
`MetaDictMap.FIT_TIER` records the tier each run was selected under, and `fit_map[key]` holds `[label, fit_dir]` with label
`Best`, `Acceptable`, `AutoQC_Best`, `AutoQC_Pass`, `Default` (or `Unreviewed` for all-fits rows).

### Data Maps
The FacetsMeta object keeps track of a variety of id -> value maps.

* master_file_dict - A map of relevant files for each run. {key: [out_file, cncf_file, qc_file, ..., ]}.  As the master file dictionary contains a large number of file types and may expand over time, please only make calls to this file dictionary using the [MetaDictMap](metadictmap.md) class.  For example, to access an adjusted .seg file, use `master_file_dict.get(id)[MetaDictMap.ADJUSTED_SEG_FILE]`
* cancer_type_map - A map of sample ids to cancer types.
* cancer_type_detail_map - A map of sample ids to cancer type detailed.
* patient_id_map - A map of sample ids to patient ids.
* clinical_purity_map - A map of purity values from the clinical sample file.
* onkotree_code_map - A map of onkotree codes.
* cvr_tmb_score_map - A map of cvr tmb scores from the clinical sample file.
* msi_score_map - A map of msi scores.
* long_id_map - A map of sample ids to their corresponding long_ids.  id -> [long_id1, long_id2...] (single-run mode: `[pair id]`).
* samples_from_file - A list of samples from a file that should be selected for this object. Populated by selectSamplesFromFile(input_file).
* fit_map - A map of run key -> [label, fit_path].
* build_counts - counters from the last build: fits per label, `multi_normal`, `missing_manifest`, `missing_class`, `bad_manifest`, `missing_files`, `no_fit_selected`, `skipped_ultra`, `skipped_not_selected`, `standard_samples`, `twon_samples`, ...
* failed_samples - run keys / pair ids whose files could not be parsed.

### Functions

* buildFacetsMeta()
  * Populates the FacetsMeta object (runs `parseClinicalSample()`), handling persistence as described above.
* parseClinicalSample()
  * Reads the clinical sample file (or walks the repository), discovers each sample's layout, reads its manifest(s) and fills master_file_dict. Raises `FacetsRepositoryError` for a missing repository path or an unreadable clinical file; per-sample problems are counted and skipped.
* setSingleRunPerSample(bool useSingleSample, bool allowDefaults=False)
  * Select a single run per sample (per class for 2n). With allowDefaults the default fit is used when no reviewed/autoQC fit exists. `False, False` selects every fit of every sample (all-fits mode).
* setFitPolicy(FitPolicy policy)
  * Fine control over single-run selection (see above).
* setVerbose(bool activateVerbose)
  * Print per-sample diagnostics (skipped fits, unreadable manifests, parse failures) during builds. Off by default.
* selectSamplesFromFile(input_file)
  * Restrict the build to the sample ids listed in a file, one per line (tumor ids `P-0000001-T01-IM6` or full pair ids; blank and `#` lines ignored).
* shortId(key) / pairIdForKey(key)
  * Id helpers, see the key table above.
* printLogo()
  * A simple method for printing a FACETS API logo. **Run automatically during buildFacetsMeta()**.

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
  * Chromosome arm position map (hg19). chrID -> [p_start, p_end, q_start, q_end]. For example `FacetsMeta.chr_arms.get(1)`.
* CN_CALL_TABLE
  * Copy-number call lookup keyed on `(wgd, tcn, lcn)`; use `FacetsRun.get_cn_call(tcn, lcn)` or `facetsAPI.get_cn_call(wgd, tcn, lcn)`.

#### Examples

```python
  # Build a FacetsMeta object, selecting one run per sample (reviewed > autoQC > default),
  # reading "purity" data, persisting the object into all_meta.dat.
  clinical_sample_file  = "/path/to/data_clinical_sample.oncokb.txt"
  facets_dir            = "/path/to/impact/facets/all/"
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "purity", "all_meta.dat")
  prepared_metadata.setSingleRunPerSample(True, True)
  prepared_metadata.buildFacetsMeta()
```

```python
  # Same, with no clinical sample file: the run is based on what is present in facets_dir.
  prepared_metadata = FacetsMeta("", facets_dir, "purity")
  prepared_metadata.setSingleRunPerSample(True, True)
  prepared_metadata.buildFacetsMeta()
```

```python
  # A 2n repository: clinical fits only, human-reviewed fits only.
  prepared_metadata = FacetsMeta(clinical_sample_file, "/path/to/impact_2n/facets/all/", "purity",
                                 fit_class="clinical", fit_policy=FitPolicy.reviewed_only())
  prepared_metadata.setSingleRunPerSample(True, False)
  prepared_metadata.buildFacetsMeta()
```

```python
  # All samples and all runs, "hisens" data, no persistence.
  prepared_metadata = FacetsMeta(clinical_sample_file, facets_dir, "hisens")
  prepared_metadata.buildFacetsMeta()
```

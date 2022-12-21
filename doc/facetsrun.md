## FacetsRun

### Initializing a FacetsRun

The FacetsRun class represents a single facets run (a.k.a. fit) as represented in a fit directory in a facets folder.
For example, a fit folder might be called /default/ or /alt_dipLogR-0.03/.
FacetsRun objects are built automatically as part of a [FacetsDataset](FacetsDataset.md), and in general should not have to be manually constructed.
A [FacetsDataset](FacetsDataset.md) will contain a set of FacetsRuns and/or FacetsSamples. A FacetsSample may contain multiple FacetsRuns.

### FacetsRun Variables
FacetsRuns contain information about the specific fit, as well as variables that may be necessary when performing
certain operations such as alteration analysis.  The following variables are available in a FacetsRun object.  

#### These values are part of a Facets Run/Fit and are extracted from files in the assocated run folder.
* segments - A list of all [FacetsSegments](FacetsSegment.md) associated with this run.
* genes - A list of all [FacetsGenes](FacetsGene.md) associated with this run.
* id - The DMP ID for this run.
* patientId - The patient ID for this run.
* fitDir - The full path to the fit directory for this run.
* cancerType - The cancer type for this run from the clinical sample data.
* cancerTypeDetail - The cancer type detail for this run from the clinical sample data.
* purity - The purity value for this run as calculated by FACETS.
* clinicalPurity - The clinical purity value as indicated in the clinical sample data.
* onkoCode - The oncoKB abbreviation code for indicating cancer type for this run as indicated in the clinical sample data.
* ploidy - The ploidy for this run.
* dipLogR - The dipLogR for this run.
* cval - the cValue for this run.
* wgd - A boolean value indicating if this sample is considered to have whole genome duplication.
* fga - The fraction of the genome that is altered in this run.
* frac_loh - The fraction of loss of heterozygosity for this run.
* facets_qc - A boolean value indicating if this sample passed facets qc.
* tmb - The tumor mutational burden for this run as indicated in the clinical sample data.
* msi - The MSI value for this run as indicated in the clinical sample data.
* purity_hisens - Whether this run was created based on purity or hisens data.  String: "purity" or "hisens".

#### These values are calculated/used in alteration analysis.
* alteration_arm_num - The number of arms that need to be LOH or Gain to make a sample level Hypoploidy or GainHeavy call.
* alteration_cov_min - The percentage of the arm that needs to be covered to make LOH/Gain calls.
* run_alterations - Whether or not to run the alteration calculation functions.
* split_arms - Whether or not to split segments into p and q arms.
* alterations - A map of {chr_arm -> [alteration, percentAltered]}, where alteration is (LoH or Gain or Neutral) based on the segments for this sample.
* isHypoploidy - A boolean value indicating if this run is considered overall to exhibit hypoploidy.
* isGainHeavy - A boolean value indicating if this run is considered overall to be gain heavy. Adjusted for WGD.
* cfLevel_arms - A map of {cf_level -> [chr_arms]}. CF levels are arm level calls that correspond to other arms with events at the same CF.
* cfLevel_values - A map of {cf_level -> cf_value}.
* arm_cf_data - A map of {arm -> [cf, percent_arm_coverage]}
* max_arm_level - The arm level that has the most members.

#### Functions
* addSegment(seg)
  * 
* removeSegment(seg)
  * 
* printAllGenes()
  *
* addGene(gene)
  *
* removeGene(gene)
  *
* getSegmentsByChromAndArm(target_chrom, target_arm)
  *
* getSegmentsByChrom(target_chrom)
  *
* calculatePurityByCF(use_base)
  *
* getNumberLossArms(percent_lost)
  *
* getNumberGainArms(percent_gain)
  *
* getSortedCFs(target_arm)
  *
* getCFbyArm()
  *
* getArmLevelCF(target_arm, match_cf_range=0.01)
  *
* mergeCFLevels(lsts)
  *
* defineCFLevels(match_cf_range=0.01)
  *
* defineArms()
  *
* defineSampleLevelAlterationCall(min_percent_arm_cov, min_arm_num)
  *
* calculateLohAndGains(selected_segs, target_chr, target_arm)
  *
* defineAllLohAndGains()
  *
* printAllSegments()
  *
* printTargetSegment(segNumber)
  *
* printSample()
  *
* parseOut(out_file)
  *
* parseGeneLevel(gene_level_file)
  *
* parseSampleQC(qc_file, expected_cval)
  *
* parseFacetsQC(sample_dir_map)
  *
* parseCNCF(cncf_file, adjseg)
  *
* buildFacetsRun(sample_id, facets_metadata)
  *

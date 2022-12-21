## FacetsRun

### Initializing a FacetsRun

The FacetsRun class represents a single facets run (a.k.a. fit) as represented in a fit directory in a facets folder.
For example, a fit folder might be called /default/ or /alt_dipLogR-0.03/.
FacetsRun objects are built automatically as part of a [FacetsDataset](facetsdataset.md), and in general should not have to be manually constructed.
A [FacetsDataset](facetsdataset.md) will contain a set of FacetsRuns and/or FacetsSamples. A FacetsSample may contain multiple FacetsRuns.

### FacetsRun Variables
FacetsRuns contain information about the specific fit, as well as variables that may be necessary when performing
certain operations such as alteration analysis.  The following variables are available in a FacetsRun object.  

#### These values are part of a Facets Run/Fit and are extracted from files in the assocated run folder.
* segments - A list of all [FacetsSegments](facetssegment.md) associated with this run.
* genes - A list of all [FacetsGenes](facetsgene.md) associated with this run.
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
  * This function will add a segment to the segment array of type [FacetsSegment](facetssegment.md). Returns True if successful.
* removeSegment(seg)
  * This function will remove a segment from the segment array of type [FacetsSegment](facetssegment.md). Returns True if successful.
* printAllGenes()
  * This function will print all [FacetsGenes](facetsgene.md) in the gene array.
* addGene(gene)
  * This function will add a [FacetsGenes](facetsgene.md) to the genes array.  Returns True if successful.
* removeGene(gene)
  * This function will remove a [FacetsGenes](facetsgene.md) from the genes array.  Returns True if successful.
* getSegmentsByChromAndArm(target_chrom, target_arm)
  * This function, given a target chromosome (int) and target arm (str), will return a list of all segments for this sample that fall in the requested regions.
* getSegmentsByChrom(target_chrom)
  * This function, given a target chromosome (int), will return a list of all segments for this sample that fall in the requested regions.
* calculatePurityByCF(use_base=False)
  * This function will calculate purity based on all of the segments associated with this sample. Purity is the maximum cf value that is not 1 or NA. If use_base is set to true (default: False), this function will calculate the purity based on the non .em cf values.  Otherwise cf.em will be used.
* getNumberLossArms(percent_lost)
  * This function will return the number of arms with LOH over 'percent_lost' in this sample.
* getNumberGainArms(percent_gain)
  * This function will return the number of arms with GAIN over 'percent_gain' in this sample.
* getSortedCFs(target_arm)
  * This function will sort the self.cfLevel_arms into and return an indexed list of arm levels. It returns two lists, 1 indexed (ignore 0 position), for [[p_arms],[q_arms]].
* getCFbyArm(target_arm)
  * This function returns cf data for a specified arm (str), i.e. "1p". 
* getArmLevelCF(target_arm, match_cf_range=0.01)
  * This function will calculate a CF value covering the largest segment of a target_arm. I.E. if an arm is 70% cf=.52 and 30% cf=.23, this function will return .52. match_cf_range is the maximum distance that two segments can be from one another in order to be considered the same CF. This function returns [longest_cf_value, percent_of_arm_covered].
* mergeCFLevels(lsts)
  * This function is for merging CF levels.  Helper for defineCFLevels(). Given a list of levels to merge [[1,3], [2,4], [2,3]] this function should merge the sub-lists into sets where any overlap between sets is given.  In this example this should return [[1,2,3,4]].
* defineCFLevels(match_cf_range=0.01)
  * This function will identify arms that belong on the same level together and assign data for CF values and maps to this FacetsRun object. Uses defined arms and determines the cf levels across all arms. A CF level is a CF value that matches other CF levels in a sample. All samples with the same CFs will be in the same level. match_cf_range is the maximum distance that two segments can be from one another in order to be called at the same level. This function will identify arms that belong on the same level together and assign data for CF values and maps to this FacetsRun object.
* defineArms()
  * Check segments to determine p or q arm. If centromere is crossed, splits the segment into two new segments. 
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

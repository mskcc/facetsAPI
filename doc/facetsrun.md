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
  * This function will check this sample and determine if it is a hypoploidy sample. It will make these determinations based on a minimum percent of arm covered for min_arm_num of samples. This function will set isHypoploidy and isGainHeavy for this FacetsRun object.
* calculateLohAndGains(selected_segs, target_chr, target_arm)
  * This is a helper function for defineAllLohAndGains().  It will perform LoH and Gain calculations on a set of segments for a given chromosome and arm.
* defineAllLohAndGains()
  * This function will process all segments for this FacetsRun object and populate the 'alterations' map.  A run is determined to be LoH if it is 0 LCN for a percentage greater than that defined in the [FacetsSegment](facetssegment.md) percent_arm_loh threshold variable, considering all segments for each arm. A gain is determined by a correspondingly TCN greater than [FacetsSegment](facetssegment.md)'s min_gain_tcn and percent_arm_gain threshold variables.  This function will process all chromosomes/arms in this FacetsRun object.
* printAllSegments()
  * Print all segments in the segment array.
* printTargetSegment(segNumber)
  * Print a specific segment by index, where segNumber is the index of the segment array to print.
* printSample()
  * Print a summary of this sample.
* parseOut(out_file)
  * This function accepts a path to a .out file from a facets directory. It will extract and return purity, cval, ploidy, and dipLogR values as a list. Returns: [purity,cval,ploidy,dipLogR]
* parseGeneLevel(gene_level_file)
  * This function accepts a path to a .gene_level.txt file and returns a list of lists containing data representing the gene level file. Returns [[gene1,start,end...],[gene2,start,end...]]
* parseSampleQC(qc_file, expected_cval)
  * This function accepts a path to a .sample.qc.txt file and an expected cval which can be identified with the parseOut function.  It returns a list containing data for whole genome duplication, fraction genome altered, and fraction LOH. Returns: [wgd,fga,frac_loh]
* parseFacetsQC(sample_dir_map)
  * This function accepts the directory map list for a sample. It will extract and return the "facets_qc" value for the target sample.
* parseCNCF(cncf_file, adjseg)
  * This function accepts a CNCF file path and parses it for data. It will build and return a list of FacetsSegment type objects representing the CNCF.
* buildFacetsRun(sample_id, facets_metadata)
  * This function accepts a sample ID and parses all relevant files and assembles all desired values.  It creates and returns an object of type FacetsRun with all manditory variables and fields populated.  It does not run alteration analysis, which can be run on all FacetsRuns in a [FacetsDataset](facetsdataset.md) by using the FacetsDataset.runAlterationAnalysis() function.

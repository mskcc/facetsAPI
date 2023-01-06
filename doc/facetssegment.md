## FacetsSegment

### Initializing a FacetsSegment

The FacetsSegment class represents a single facets segment as represented in a facets run's .seg files. 
FacetsSegment objects are built automatically as part of a [FacetsRun](facetsrun.md), and in general should not have to be manually constructed.  

### FacetsSegment Variables
The following variables are available in a FacetsSegment object.  

#### These values are used in or calculated in FacetsRun.defineArms().
* percent_arm_loh - The percentage of the arm that needs to be LCN=0 to make an arm level LOH call.
* min_gain_tcn - The minimum value of TCN required to consider a segment to be a Gain call.
* percent_arm_gain - The percentage of the arm that needs to be TCN=min_gain_tcn to make an arm level Gain call.
* arm - The chr/arm.  I.E. 1p, 5q, etc.
* percentArm - The percentage of the arm that this segment covers. 
* isLoH - Boolean. If this segment is considered a LoH segment
* isGain - Boolean. If this segment is considered a Gain segment.

#### These values are extracted from the .CNCF file unless otherwise indicated.
* chrom - Chromosome number.
* star - Segment start position.
* end - Segment end position
* cnlr_median - Unadjusted mean values.  These correspond to the mean column in .unadjusted.seg files.
* length - length of this segment.
* cf - Clonal fraction for this segment.  Note that this is the expectation maximization adjusted cf value (em).
* cf_base - The base clonal fraction.  Unadjusted.
* num_mark - The number of SNPs in the segment
* nhet - Minimum number of heterozygote snps in a segment used to call minor cn
* segclust - The segment cluster to which segment belongs.
* cnlr_median_clust - The median log-ratio of the segment.
* mafR - The log-odds-ratio summary for the segment
* adj_mean - Adjusted mean values.  These correspond to the mean column in .adjusted.seg files.
* mafR_clust - The log-odds-ratio summary for the segment cluster.
* tcn - Total copy number for this segment. Note that this is the expectation maximization adjusted cf value (em).
* lcn - Lesser copy number for this segment. Note that this is the expectation maximization adjusted cf value (em).


#### Functions
* compareSegments(segToCompare):
  * This function will compare a segment with the segment provided as a parameter.  If they have the same start, end, and chrom, this function returns True.  
* printSegment()
  * This function will print a report for the given segment. 

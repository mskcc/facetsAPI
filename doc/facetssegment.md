## FacetsSegment

### Initializing a FacetsSegment

The FacetsSegment class represents a single facets segment as represented in a facets run's .seg files. 
FacetsSegment objects are built automatically as part of a [FacetsRun](FacetsRun.md), and in general should not have to be manually constructed.  

### FacetsSegment Variables
The following variables are available in a FacetsSegment object.  

#### These values are used in or calculated in FacetsRun.defineArms().
* percent_arm_loh - The percentage of the arm that needs to be LCN=0 to make an arm level LOH call.
* min_gain_tcn - The minimum value of TCN required to consider a segment to be a Gain call.
* percent_arm_gain - The percentage of the arm that needs to be TCN=min_gain_tcn to make an arm level Gain call.
* self.arm - The chr/arm.  I.E. 1p, 5q, etc.
* self.percentArm - The percentage of the arm that this segment covers. 
* self.isLoH - Boolean. If this segment is considered a LoH segment
* self.isGain - Boolean. If this segment is considered a Gain segment.

#### These values are extracted from the .CNCF file unless otherwise indicated.
* self.chrom - Chromosome number.
* self.star - Segment start position.
* self.end - Segment end position
* self.cnlr_median - Unadjusted mean values.  These correspond to the mean column in .unadjusted.seg files.
* self.length - length of this segment.
* self.cf - Clonal fraction for this segment.  Note that this is the expectation maximization adjusted cf value (em).
* self.cf_base - The base clonal fraction.  Unadjusted.
* self.num_mark - **ALLISON**
* self.nhet - **ALLISON**
* self.segclust - **ALLISON**
* self.cnlr_median_clust - **ALLISON**
* self.mafR - **ALLISON**
* self.adj_mean - Adjusted mean values.  These correspond to the mean column in .unadjusted.seg files.
* mafR_clust - **ALLISON**
* tcn - Total copy number for this segment. Note that this is the expectation maximization adjusted cf value (em).
* lcn - Lesser copy number for this segment. Note that this is the expectation maximization adjusted cf value (em).


#### Functions

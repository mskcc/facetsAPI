## FacetsGene

### Initializing a FacetsGene

The FacetsGene class represents a single facets gene as represented in a facets run's gene_level.txt file.  
FacetsGene objects are built automatically as part of a [FacetsRun](FacetsRun.md), and in general should not have to be manually constructed.  

#### FacetsGene Variables
The following variables are available in a FacetsGene object.  

* gene - The gene name for this FacetsGene object.
* gene_start - The start position for this gene.
* gene_end - The end position for this gene.
* seg_start - The start position of the FACETS segment associated with this gene.
* seg_end - The end position of the FACETS segment associated with this gene.
* seg_length - The length of the FACETS segment associated with this gene.
* cf - The clonal fraction of the FACETS segment associated with this gene.  Note this value reflects cf.em, not the uncorrected base cf value.
* tcn - The total copy number of the FACETS segment associated with this gene. tcn.em value.
* lcn - The lesser copy number of the FACETS segment associated with this gene. lcn.em value.
* cn_state - The copy number state of the FACETS segment associated with this gene as called by FACETS.
* filter - A boolean value indicating if this segment passed filter criteria.  



#### Functions
* printGene()
  * This function will print a report of the details of this gene.
  
* compareGenes(geneToCompare)
  * This function will compare two genes to see if they are the same.  Comparison is based on gene name, gene start, and gene end position.  Returns true if the genes match. 


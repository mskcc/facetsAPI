
library(data.table)
library(doParallel)
library(facetsPreview, lib.loc = '/juno/work/ccs/pricea2/ondemand/R_libs/')

#options(echo=TRUE) # To see commands in output file

args <- commandArgs(trailingOnly = TRUE)

facets_preview_config_file = args[1]
tumor_sample = args[2]
sample_id = args[3]
sample_path = args[4]
cbio_maf_path = args[5]
cbio_maf_nonsignedout_path = args[6]

cbio_maf <- fread(cbio_maf_path)
cbio_maf_nonsignedout <- fread(cbio_maf_nonsignedout_path)

#tumor_sample
#Tumor_Sample_Barcode
#P-0000280-T05-IM7

#sample_id
#P-0071250-T01-IM7_P-0071250-N01-IM7

#print(facets_preview_config_file)
#print(tumor_sample)
#print(sample_id)
#print(sample_path)
#print(cbio_maf_path)
#print(cbio_maf_nonsignedout_path)

write.table(
  cbio_maf %>% 
    filter(Tumor_Sample_Barcode == tumor_sample),
  file=paste0(sample_path, '/', sample_id, '.maf'),
  quote=F, row.names=F, sep='\t')

 system(paste0('chmod -R 775 ', sample_path, '/', sample_id, '.maf'))

write.table(
  cbio_maf_nonsignedout %>% 
    filter(Tumor_Sample_Barcode == tumor_sample),
  file=paste0(sample_path, '/', sample_id, '.nonsignedout.maf'),
  quote=F, row.names=F, sep='\t')

  system(paste0('chmod -R 775 ', sample_path, '/', sample_id, '.nonsignedout.maf'))


facetsPreview::generate_genomic_annotations(sample_id, sample_path, facets_preview_config_file, T)

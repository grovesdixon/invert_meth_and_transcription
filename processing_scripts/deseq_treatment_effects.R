#deseq_treatment_effects.R 
library(DESeq2)

#SET UP THE DATA TO RUN DESEQ
rm(list=ls())
source('invert_meth_ge_functions.R')
bioprojectDirs




# LOAD THE DESEQ INPUT FILES FOR INVERTS ----------------------------------------------

deseq_path_list = paste('bioprojects/', bioprojectDirs, '/rna_results/deseqInput.Rdata', sep='')
names(deseq_path_list) = dat_names


run_deseq = function(n){
  print(paste('Running DEseq on ', n, '...', sep=''))
  ll = load(deseq_path_list[[n]])
  dds<-DESeqDataSetFromMatrix(counts,
                              colData = coldata, 
                              design = formula(~ treat))
  dds <- DESeq(dds)
  res = results(dds, contrast = c('treat', 'treated', 'control'))
  return(res)
}

res_list = map(dat_names, run_deseq)
names(res_list) = dat_names
save(res_list, file='figure_plotting/deseq_res_list.Rdata')


# REPEAT FOR PRIMATES ----------------------------------------------
source('metadata/primate_bioproject_source.R')
bioprojectDirs
dat_names
deseq_path_list = paste('bioprojects/', bioprojectDirs, '/rna_results/deseqInput.Rdata', sep='') #build these with initialize_feature_counts.R in each species' subdirectory
names(deseq_path_list) = dat_names


run_deseq_liver = function(n){
  print(paste('Running DEseq on ', n, '...', sep=''))
  ll = load(deseq_path_list[[n]])
  coldata$liver = if_else(coldata$tissue=='liver',
                          'yes',
                          'no')
  dds<-DESeqDataSetFromMatrix(counts,
                              colData = coldata, 
                              design = formula(~ liver + Individual))
  dds <- DESeq(dds)
  res = results(dds, contrast = c('liver', 'yes', 'no'))
  return(res)
}

liver_list = map(dat_names, run_deseq_liver)
names(liver_list) = dat_names
save(liver_list, file='figure_plotting/primate_deseq_res_list.Rdata')

#deseq_treatment_effects.R 

#SET UP THE DATA TO RUN DESEQ
rm(list=ls())
source('invert_meth_ge_functions.R')
bioprojectDirs




# LOAD THE DESEQ INPUT FILES ----------------------------------------------

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

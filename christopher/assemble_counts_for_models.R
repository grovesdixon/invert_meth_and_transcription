#assemble_counts_for_models.R

rm(list=ls())
source('invert_meth_ge_functions.R')
source('metadata/species_logos/logo_source.R')

#add amillepora to bioprojectDirs
bioprojectDirs['amillepora'] = 'amillepora_PRJNA601565'
dat_names = append(dat_names, 'amillepora')

#GLOBAL VARS
bioprojectDirs
dat_names



# LOAD THE GBM TOTALS AND POSITIONS --------------------------------------------------

#set up the path list
gbm_path_list = paste('bioprojects/', bioprojectDirs, '/meth_response/gene_change_methylKit.Rdata', sep='')
names(gbm_path_list) = dat_names

#load the summed gbm counts and position info, filter away genes with no counts
load_gbm_totals_and_pos = function(file_path){
  print(file_path)
  ll=load(file_path)
  d = wbounds %>% 
    left_join(totCounts, by = c("chr", "start", "end", "name")) %>% 
    group_by(chr, start, end, name) %>% 
    summarize(numMeth=sum(nC, na.rm=TRUE),
              numUnmeth=sum(nT, na.rm=TRUE)) %>% 
    filter((numMeth + numUnmeth) > 0) %>% 
    as_tibble()
}

#run
gbm_list = map(gbm_path_list, load_gbm_totals_and_pos)
names(gbm_list)
map(gbm_list, head)

#get just the position info
pos_list = map(gbm_list, function(x) return(x %>% 
                                              dplyr::select(chr:name)))

# LOAD THE GBM COUNTS ---------------------------------------------------

#set up the path list
gene_path_list = paste('bioprojects/', bioprojectDirs, '/meth_region_counts/gene_change_methylKit.Rdata_methylKit_regionCounts.Rdata', sep='')
names(gene_path_list) = dat_names

#funciton to read in rld data
load_bs_counts = function(file_path){
  print(paste('loading ', file_path, '...', sep=''))
  ll=load(file_path)
  return(my.reg.counts)
}

mkit_counts = map(gene_path_list, load_bs_counts)
map(mkit_counts, head)

# LOAD TREATMENT TABLES ---------------------------------------------------

mtt_path_list = paste('bioprojects/', bioprojectDirs, '/bisulfite_treat_table.txt', sep='')
names(mtt_path_list) = dat_names
mtt_list = map(mtt_path_list, read_tsv)
map(mtt_list, head)

#doublecheck the treatment tables line up with methylkit results files
for (n in dat_names){
  check = sum(attr(mkit_counts[[n]], 'treatment')==mtt_list[[n]]$treat)==nrow(mtt_list[[n]])
  print(paste(n, 'lines up =', check))
}


# FORMAT THE COUNTS WITH TREATMENT INFO -----------------------------------

format_methylkit_counts = function(n){
  #format the counts
  mkit = mkit_counts[[n]]
  pos_df = pos_list[[n]]
  
  #pull out data and merge with positions
  mkit_df = mkit %>% 
    data.frame() %>% 
    dplyr::select(-starts_with('coverage'),
                  -strand) %>% 
    right_join(pos_df, by=c('chr', 'start', 'end')) %>% 
    as_tibble()
  
  #parse the treatment data from attribute of the methylkit object
  treats = attr(mkit, 'treatment')
  tdf = tibble(num = 1:length(treats),
               treat = treats)
  treated_nums = tdf %>% 
    filter(treat==1) %>% 
    pull(num)
  control_nums = tdf %>% 
    filter(treat==0) %>% 
    pull(num)
  
  #isolate the methylated and unmethylated counts and pivot longer
  pos_cols = c('chr', 'start', 'end', 'name')
  meth_counts = mkit_df %>% 
    dplyr::select(pos_cols,
                  contains('numCs')) %>% 
    pivot_longer(contains('numCs'),
                 names_to = 'sample_num',
                 values_to = 'numMeth',
                 names_prefix = 'numCs')
  
  unmeth_counts = mkit_df %>% 
    dplyr::select(pos_cols,
                  contains('numTs')) %>% 
    pivot_longer(contains('numTs'),
                 names_to = 'sample_num',
                 values_to = 'numUnmeth',
                 names_prefix = 'numTs')
  
  #merge back up
  counts_df = meth_counts %>% 
    inner_join(unmeth_counts, by = c(pos_cols, 'sample_num')) %>% 
    mutate(treatment = if_else(sample_num %in% as.character(treated_nums),
                               '1',
                               '0'))
}

#run formatting
meth_data_list = list()
for (n in dat_names){
  print(n)
  mc = format_methylkit_counts(n)
  meth_data_list[[n]] = mc
}
names(meth_data_list)
map(meth_data_list, head)


#check bimodal
bm_list = list()
for (n in dat_names){
  mc = meth_data_list[[n]]
  plt = mc %>% 
    group_by(name) %>% 
    summarize(totMeth = sum(numMeth, na.rm=TRUE),
              totUnmeth = sum(numUnmeth, na.rm=TRUE),
              l.fracMeth = log(totMeth/(totMeth+totUnmeth), 2)) %>% 
    ggplot(aes(x=l.fracMeth)) +
    geom_histogram() +
    labs(subtitle=n)
  bm_list[[n]] = plt
}
plot_grid(plotlist = bm_list, nrow=3)


# LOAD DESEQ RESULTS  ---------------------------------------------------
ll=load('figure_plotting/deseq_res_list.Rdata')
names(res_list)

#add amil ge results
ll=load('bioprojects/amillepora_PRJNA601565/rna_results/tissue_results.Rdata')
res_list[['amillepora']] = res

#function to format DESeq results object
pull_res_data = function(r){
  r %>% 
    data.frame() %>% 
    rownames_to_column('name') %>% 
    dplyr::select(name, log2FoldChange, pvalue) %>% 
    as_tibble()
}

ge_data_list = map(res_list, pull_res_data)
names(ge_data_list)
map(ge_data_list, head)



# LOAD THE 1KB WINDOW DATA ------------------------------------------------

#set up the path list
window_path_list = paste('bioprojects/', bioprojectDirs, '/meth_response/significant/1KbWindows_change_significant_closeGeneDistances.tsv', sep='')
names(window_path_list) = dat_names
window_path_list = window_path_list[dat_names != 'subcaste']

#skip subcaste because has no significant windows

#load the summed gbm counts and position info, filter away genes with no counts
read_window = function(file_path){
  print(file_path)
  wdat = read_tsv(file_path)
  description = sub('ID=', '', wdat$gene.description)
  gene_name = sapply(description, function(x) strsplit(x,';')[[1]][1])
  wdat %>% 
    mutate(name = gene_name) %>% 
    dplyr::select(-gene.description, -strand)
}

window_data_list = map(window_path_list, read_window)
names(window_data_list) = dat_names[dat_names != 'subcaste']
map(window_data_list, head)


# SAVE EVERYTHING ---------------------------------------------------------

save(mtt_list, file='christopher/mtt_path_list.Rdata')
save(meth_data_list, file='christopher/meth_data_list.Rdata')
save(ge_data_list, file = 'christopher/ge_data_list.Rdata')
save(window_data_list, file = 'christopher/window_data_list.Rdata')






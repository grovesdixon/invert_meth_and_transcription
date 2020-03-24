#plot_pcas.R

rm(list=ls())
library(DESeq2)
library(methylKit)
source('invert_meth_ge_functions.R')
source('metadata/species_logos/logo_source.R')


#GLOBAL VARS
bioprojectDirs
NTOP = 10000
FIX_COORDS = FALSE


# LOAD THE RNA RLD DATA ---------------------------------------------------

#set up the path list
pathList = paste('bioprojects/', bioprojectDirs, '/rna_results/rld.Rdata', sep='')
names(pathList) = dat_names

#funciton to read in rld data
load_rld = function(x){
  print(paste('loading ', x, '...', sep=''))
  ll=load(x)
  return(list('rld'=rld.df,
              'coldata'=coldata))
}

#get the rld data for each project
rld_list = map(pathList, load_rld)

#funnction to plot the PCA
rna_pca_list = list()
for (n in dat_names){
  x=rld_list[[n]]
  rld = x[['rld']]
  coldata = x[['coldata']]
  pca_df = build_pca(rld, coldata, ntop = NTOP, pcs = 2)
  plt = plot_rld_pca(pca_df,
                     group_col = 'treatInfo',
                     fix_coords = FIX_COORDS) + 
    theme(legend.position = 'right') +
    labs(subtitle = 'RNA')
  rna_pca_list[[n]] = plt
}

plot_grid(plotlist = rna_pca_list)



# PLOT GBM PCA ------------------------------------------------------------

#set up the path list
rc_paths = paste('bioprojects/', bioprojectDirs, '/meth_region_counts/gene_change_methylKit.Rdata_methylKit_regionCounts.Rdata', sep='')
tt_paths = paste('bioprojects/', bioprojectDirs, '/bisulfite_treat_table.txt', sep='')
names(rc_paths) = names(bioprojectDirs)
names(tt_paths) = names(bioprojectDirs)

#upload the gbm proportions and treatment tables
m_prop_list = list()
tt_list = list()
for (n in names(bioprojectDirs)){
  print(paste(n,'...',sep=''))
  rc_file = rc_paths[[n]]
  tt_file = tt_paths[[n]]
  mdat = get_meth_proportion_from_reg_counts(rc_file, tt_file)
  m_prop_list[[n]] = mdat
  tt_list[[n]] = read_tsv(tt_file)
}

#BUILD THE PCA DF LIST
pca_df_list = list()
for (n in names(bioprojectDirs)){
  print(n)
  tt=tt_list[[n]]
  m_df = m_prop_list[[n]] %>% 
    tidyr::unite('reg', chr,start,end, sep='_') %>% 
    filter(!duplicated(reg)) %>% 
    column_to_rownames('reg') %>% 
    na.omit()
  pca_df = build_pca(m_df, tt, ntop = NTOP, pcs = 2)
    pca_df_list[[n]] = pca_df
}

#BUILD THE PCA PLOT LIST
meth_pca_list = list()
for (n in names(bioprojectDirs)){
  print(n)
  plt = plot_rld_pca(pca_df_list[[n]],
               group_col = 'treatInfo',
               fix_coords = FIX_COORDS) + 
    theme(legend.position = 'right') +
    labs(subtitle='GBM')
  meth_pca_list[[n]] = plt
  
}

plot_grid(plotlist = meth_pca_list)




# BUILD FINAL PCAS TOGETHER -----------------------------------------------

#some plotting vars
LOGO_CORNER = 0.95
LOGO_SCALE = 0.5
BUFFER_LEN = 0.125


plot_paired_pca = function(n){
  logo = logo_list[[n]]
  legend = cowplot::get_legend(rna_pca_list[[n]] )
  plt1 = rna_pca_list[[n]] + theme(legend.position = 'none')
  plt2 = meth_pca_list[[n]] + theme(legend.position = 'none')
  lplt = ggdraw(legend) + draw_image(logo,
                                      x=0,
                                      y=0.5,
                                      width = LOGO_SCALE,
                                      height = LOGO_SCALE)
  
  plot_grid(plt1, plt2, lplt, nrow=1, rel_widths = c(1,1,0.6))
}

paired_plt_list = map(dat_names, plot_paired_pca)
length(paired_plt_list)
plot_grid(plotlist = paired_plt_list[1:4], nrow=4)
plot_grid(plotlist = paired_plt_list[5:8], nrow=4)





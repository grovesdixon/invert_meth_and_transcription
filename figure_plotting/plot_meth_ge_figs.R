#plot_meth_level_figs.R
library(ggimage)
library(ggtree)

rm(list=ls())
source('invert_meth_ge_functions.R')
source('metadata/species_logos/logo_source.R')


#GLOBAL VARS
bioprojectDirs
SET_BREAKS = log(c(2^SET_ZERO, 0.02, 0.15, 1), 2) #X axis break points for each panel


# LOAD THE GBM LEVEL DATA ---------------------------------------------------

#set up the path list
gene_path_list = paste('bioprojects/', bioprojectDirs, '/meth_level/gene_basicStatsBed.tsv', sep='')
exon_path_list = paste('bioprojects/', bioprojectDirs, '/meth_level/exon_basicStatsBed.tsv', sep='')
promoter_path_list = paste('bioprojects/', bioprojectDirs, '/meth_level/promoter_basicStatsBed.tsv', sep='')
names(gene_path_list) = dat_names
names(exon_path_list) = dat_names
names(promoter_path_list) = dat_names

#funciton to read in rld data (note setting zero at SET_ZERO)
load_lvl = function(x){
  print(paste('loading ', x, '...', sep=''))
  df1 = read_tsv(x)
  minNonZero = df1 %>% 
    filter(fracMeth > 0) %>% 
    pull(fracMeth) %>% 
    min()
  df1 %>% 
    mutate(fracMeth = if_else(fracMeth ==0,
                              minNonZero,
                              fracMeth)) %>% 
    group_by(chr, name) %>% 
    mutate(lfrac_meth = log(fracMeth, 2),
           lfrac_meth = if_else(lfrac_meth < SET_ZERO,
                                       SET_ZERO,
                                lfrac_meth)) %>% 
    dplyr::select(-fracMeth)
}

#get the rld data for each project
gene_lvl_list = map(gene_path_list, load_lvl)
exon_lvl_list0 = map(exon_path_list, load_lvl)
promoter_lvl_list = map(promoter_path_list, load_lvl)

#SWAP OUT GENE NAMES FOR EXON LEVEL DFS AND TAKE AVERAGE

swap_exon_names = function(n){
  print(n)
  bpd = bioprojectDirs[n]
  cds_gene_file = paste('bioprojects/', bpd, '/cdsID_geneID.tsv', sep='')
  cds_gene = read_tsv(cds_gene_file) %>% 
    dplyr::rename(name=cds)
  exon_df = exon_lvl_list0[[n]]
  exon_df %>% 
    left_join(cds_gene, by = 'name') %>% 
    ungroup() %>% 
    mutate(name=gene) %>% 
    dplyr::select(-gene) %>% 
    group_by(name) %>% 
    summarize(lfrac_meth = mean(lfrac_meth))
}

exon_lvl_list = lapply(dat_names, function(x) swap_exon_names(x))
names(exon_lvl_list) = dat_names

# PLOT HISTOGRAMS ---------------------------------------------------------

#LOAD THE LOGOS
source('metadata/species_logos/logo_source.R')
logo_corner = 0.95
logo_scale = 0.2

#BUILD HISTOGRAMS
gene_hist_list = list()
exon_hist_list = list()



for (n in dat_names){
  print(n)
  logo = logo_list[[n]]
  gdf = gene_lvl_list[[n]]
  edf = exon_lvl_list[[n]]
  ghist = plot_bimodal_logo(gdf, logo, logo_scale, BREAKS=SET_BREAKS)
  ehist = plot_bimodal_logo(edf, logo, logo_scale)
  gene_hist_list[[n]] = ghist
  exon_hist_list[[n]] = ehist
}

#PLOT

assemble_row_panels = function(plot_list, nrow, xlab, ylab, relYlab, relXlab){
  pxlab = ggdraw() + draw_label(xlab)
  pylab = ggdraw() + draw_label(ylab, angle=90)
  pans = plot_grid(plotlist = plot_list, nrow=nrow)
  top = plot_grid(pylab, pans, nrow=1, rel_widths = c(relYlab, 1))
  hists = plot_grid(top, pxlab, nrow=2, rel_heights = c(1,relXlab))
}

ry = 1/30
rx = 1/15
ghists = assemble_row_panels(gene_hist_list, 2, '%methylation', 'gene count', ry, rx)
ehists = assemble_row_panels(exon_hist_list, 2, '%methylation', 'gene count', ry, rx)


# BUILD CLADOGRAM TREE ----------------------------------------------------------

#BUILD TREE
#read in the cladogram
nwk <- 'metadata/species_logos/cladogram.newick'
tree <- read.tree(nwk)
tplt = ggtree(tree, lwd=1)

#BUILD LOGOS FOR TREE
tree_logo_scale = 0.85
tree_logo_list = list()
for (n in dat_names){
  print(n)
  logo = logo_list[[n]]
  lplt = ggdraw() + 
    draw_image(logo,
               x=0,
               y=.1,
               width = tree_logo_scale,
               height = tree_logo_scale)
  tree_logo_list[[n]] = lplt
}

#BUILD HISTOGRAMS
tree_hist_list = list()
for (n in dat_names){
  print(n)
  gdf = gene_lvl_list[[n]]
  thist = plot_bimodal(gdf, BREAKS=SET_BREAKS) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank())
  tree_hist_list[[n]] = thist
}


# DESCRIPTION PANELS ------------------------------------------------------

#set up the path list
library(readxl)
ddat = read_excel('bioprojects/project_features.xlsx')

#build description list
desc_list = list()
for (d in bioprojectDirs){
  print(d)
  sub_ddat = ddat %>% 
    filter(bioproject==d)
  print(sub_ddat)
  strings = as.character(sub_ddat$value)
  names(strings) = as.character(sub_ddat$feature)
  strings = strings[c('species', 'reference', 'bioproject', 'treatment', 'methN', 'rnaN')]
  single = paste(strings, collapse = '\n')
  strings_plt = ggdraw() + draw_label(single, hjust=0, x=0.1, size=11)
  desc_list[[d]] = strings_plt
}


# ASSEMBLE THE TREE -------------------------------------------------------

#tree part
logo_panels = plot_grid(plotlist = rev(tree_logo_list), nrow=length(tree_logo_list))
tree_plt = plot_grid(tplt,logo_panels)
hist_panels = plot_grid(plotlist = rev(tree_hist_list), nrow=length(gene_hist_list))
desc_panels = plot_grid(plotlist = rev(desc_list), nrow=length(desc_list))
plot_grid(tree_plt,
          hist_panels,
          desc_panels,
          nrow=1, rel_widths = c(1.3,1,1.5))


# GET FPKM ----------------------------------------------------------------

#GET GENE LENGTHS

#function to pull gene lengths
get_gene_lengths = function(x){
  print(paste('loading ', x, '...', sep=''))
  df = read_tsv(x) %>% 
    mutate(length = end - start) %>% 
    dplyr::select(name, length)
  return(df)
}

gene_length_list = map(gene_path_list, get_gene_lengths)
names(gene_length_list) = dat_names

#GET FPKM

#set up the path list
counts_path_list = paste('bioprojects/', bioprojectDirs, '/rna_results/deseqInput.Rdata', sep='')
names(counts_path_list) = dat_names

#function to get fpkm from captured fraction counts
get_fpkm = function(counts, length_df){
  dat = counts %>% 
    rownames_to_column('name') %>% 
    left_join(length_df,
              by='name')
  lengths = dat$length
  m = apply(counts, 2, function(x) sum(x)/1e6)
  k=lengths/1e3
  fpkm=sweep(counts, 1, k, `/`) %>% 
    sweep(2, m, `/`)
  return(fpkm)
}

#funciton to read in rld data
get_fpkm_for_dat = function(n){
  print(n)
  counts_path = counts_path_list[[n]]
  length_df = gene_length_list[[n]]
  ll=load(counts_path)
  fpkm_df = get_fpkm(counts, length_df)
  mn_fpkm = apply(fpkm_df, 1, function(x) mean(x, na.rm=TRUE))
  res_df = tibble('name' = rownames(fpkm_df),
               'lfpkm' = log(mn_fpkm, 2))
  return(res_df)
}

#get mean fpkm for each dataset
mnfpkm_list = list()
for (n in dat_names){
  mnfpkm_list[[n]] = get_fpkm_for_dat(n)
}



# PLOT METH LVL AND GE LVL ------------------------------------------------

#function to build the plot
plot_lvl_lvl = function(n, mlvl_list){
  logo = logo_list[[n]]
  m_df = mlvl_list[[n]]
  r_df = mnfpkm_list[[n]]
  merge_df = inner_join(m_df, r_df,
                        by='name')
  plt = plot_scatter_pearsonCor_annotated(merge_df, 'lfrac_meth', 'lfpkm', 'GBM', 'expression level', ALPHA=0.1) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(breaks = SET_BREAKS,
                       labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank())
    ggdraw(plt) + draw_image(logo,
                              x=logo_corner - logo_scale,
                              y=logo_corner-logo_scale, 
                              width = logo_scale,
                              height = logo_scale)
}

#run accross the datasets
glvls_plt_list = list()
elvls_plt_list = list()
plvls_plt_list = list()
for (n in dat_names){
  glvls_plt_list[[n]] = plot_lvl_lvl(n, gene_lvl_list)
  elvls_plt_list[[n]] = plot_lvl_lvl(n, exon_lvl_list)
  plvls_plt_list[[n]] = plot_lvl_lvl(n, promoter_lvl_list)
}


#assmeble multiplanels
ry = 1/30
rx = 1/15
glvl_ge = assemble_row_panels(glvls_plt_list, 2, 'gbM % methylation', 'expression level', ry, rx)
elvl_ge = assemble_row_panels(elvls_plt_list, 2, 'exon % methylation', 'expression level', ry, rx)
plvl_ge = assemble_row_panels(elvls_plt_list, 2, 'promoter % methylation', 'expression level', ry, rx)




# GBM LVL VS PROMOTER METHYLATION LVL -------------------------------------

#promter methylatio shows almost exact same pattern as GBM
#so plot them against each other to show they are basically the same

gene_lvl_list = map(gene_path_list, load_lvl)
promoter_lvl_list = map(promoter_path_list, load_lvl)

head(gene_lvl_list[[n]])
head(promoter_lvl_list[[n]])


#assemble a list of merged gbm and promoter levels
gbm_promoter_lvls = list()
for (n in dat_names){
  print(n)
  pdat = promoter_lvl_list[[n]] %>% 
    dplyr::rename(promoter_lvl = lfrac_meth)
  gdat = gene_lvl_list[[n]] %>% 
    dplyr::rename(gbm_lvl = lfrac_meth)
  merge_df = pdat %>% 
    left_join(gdat, by = c('chr', 'name'))
  gbm_promoter_lvls[[n]] = merge_df
}


#plot their relationship
gbm_promoter_plts=list()
for (n in dat_names){
  merge_df = gbm_promoter_lvls[[n]]
  logo = logo_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(merge_df,
                                          'gbm_lvl','promoter_lvl', 
                                          'gbM % methylation', 'promoter % methylation', 
                                          ALPHA=0.1, ylim=c(SET_ZERO,0.1)) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(labels = log2_to_percent) +
    scale_y_continuous(labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank()) +
    geom_vline(xintercept = log(0.025, 2), color='red', lty=2) +
    geom_hline(yintercept = log(0.025, 2), color='red', lty=2)
  plt2 = ggdraw(plt) + draw_image(logo,
                                  x=logo_corner - logo_scale*3,
                                  y=logo_corner-logo_scale, 
                                  width = logo_scale,
                                  height = logo_scale)
  gbm_promoter_plts[[n]] = plt2
}


#assmeble multiplanels
ry = 1/30
rx = 1/15
glvl_plvl = assemble_row_panels(gbm_promoter_plts, 2, 'gbM % methylation', 'promoter % methylation', ry, rx)
glvl_plvl


# PROMOTER LEVEL VS GE LEVEL ----------------------------------------------
map(gbm_promoter_lvls, head)



# PROMOTER LEVEL AND GE LEVEL BY GBM-PROMOTER GROUPING --------------------

#first assign group types:
# M body M promoter
# M body Um promoter
# Um body M promoter
# Um body Um promoter

head(gbm_promoter_lvls[[n]])

assign_epistates = function(df){
  lcut = log(cut, 2)
  df %>% 
    dplyr::select(chr, name, gbm_lvl, promoter_lvl) %>% 
    mutate(Mg = gbm_lvl > lcut,
           Mp = promoter_lvl > lcut,
           epi_class = 'unassigned',
           epi_class = if_else(Mg & Mp,
                               'Mg.Mp',
                               epi_class),
           epi_class = if_else(Mg & !Mp,
                               'Mg.Up',
                               epi_class),
           epi_class = if_else(!Mg & Mp,
                               'Ug.Mp',
                               epi_class),
           epi_class = if_else(!Mg & !Mp,
                               'Ug.Up',
                               epi_class),
           gbm_class = if_else(Mg,
                               'methylated',
                               'unmethylated'))
}
cut= 0.025
epi_state_list = map(gbm_promoter_lvls, assign_epistates)
names(epi_state_list)

promoter_gelvl_by_gbm = list()
for (n in dat_names){
  e_df = epi_state_list[[n]]
  r_df = mnfpkm_list[[n]]
  logo = logo_list[[n]]
  merge_df = inner_join(e_df, r_df,
                        by='name')
  
  plt = merge_df %>% 
    filter(!is.na(gbm_class)) %>% 
    ggplot(aes(x=promoter_lvl, y=lfpkm)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method='lm') +
    scale_x_continuous(labels = log2_to_percent) +
    facet_wrap(~epi_class)
  plt2 = ggdraw(plt) + draw_image(logo,
                           x=logo_corner - logo_scale,
                           y=logo_corner-logo_scale, 
                           width = logo_scale,
                           height = logo_scale)
  promoter_gelvl_by_gbm[[n]] = plt2
}

plot_grid(plotlist = promoter_gelvl_by_gbm, nrow=2)


# PROMOTER LEVEL AND GE LEVEL AMONG METHYLATED GENES ----------------------

#plot their relationship
Mgbm_promoter_plts=list()
for (n in dat_names){
  e_df = epi_state_list[[n]]
  r_df = mnfpkm_list[[n]]
  logo = logo_list[[n]]
  merge_df = inner_join(e_df, r_df,
                        by='name') %>% 
    filter(gbm_class=='methylated') 
  plt = plot_scatter_pearsonCor_annotated(merge_df,
                                          'promoter_lvl','lfpkm', 
                                          'promoter % methylation', 'GE level', 
                                          ALPHA=0.1) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank())
  plt2 = ggdraw(plt) + draw_image(logo,
                                  x=logo_corner - logo_scale*3,
                                  y=logo_corner-logo_scale, 
                                  width = logo_scale,
                                  height = logo_scale)
  Mgbm_promoter_plts[[n]] = plt2
}


#assmeble multiplanels
ry = 1/30
rx = 1/15
Mglvl_plvl = assemble_row_panels(Mgbm_promoter_plts, 2, 'promoter % methylation', 'mRNA level', ry, rx)
Mglvl_plvl



# PLOT DIFFERENTIAL EXPRESSION BETWEEN GROUPS -----------------------------

ll=load('figure_plotting/deseq_res_list.Rdata') #use processing_scripts/deseq_treatment_effects.R to build this
ll

plot_lvl_change = function(n, mlvl_list){
  logo = logo_list[[n]]
  m_df = mlvl_list[[n]]
  r_df = res_list[[n]] %>% 
    data.frame() %>% 
    rownames_to_column('name') %>% 
    mutate(abs_diff = abs(log2FoldChange))
  merge_df = inner_join(m_df, r_df,
                        by='name')
  plt = plot_scatter_pearsonCor_annotated(merge_df, 'lfrac_meth', 'abs_diff', 'GBM', 'expression difference',
                                          ALPHA=0.1,
                                          ylim = c(0,4)) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(breaks=SET_BREAKS,
                       labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(limits=c(0,4)) +
    theme(axis.title = element_blank())
  ggdraw(plt) + draw_image(logo,
                           x=logo_corner - logo_scale,
                           y=logo_corner-logo_scale, 
                           width = logo_scale,
                           height = logo_scale)
}


#run the function
g_change_plt_list = list()
e_change_plt_list = list()
for (n in dat_names){
  g_change_plt_list[[n]] = plot_lvl_change(n, gene_lvl_list)
  e_change_plt_list[[n]] = plot_lvl_change(n, exon_lvl_list)
  
}


#assemble multiplanel
ry = 1/30
rx = 1/15
gdiff_ge = assemble_row_panels(g_change_plt_list, 2, '% methylation', 'expression difference', ry, rx)
ediff_ge = assemble_row_panels(e_change_plt_list, 2, '% methylation', 'expression difference', ry, rx)


# ASSEMBLE GE FIGURES -----------------------------------------------

#other panels
lvl_panels =  plot_grid(plotlist = glvls_plt_list, nrow=length(glvls_plt_list))
diff_panels =  plot_grid(plotlist = g_change_plt_list, nrow=length(g_change_plt_list))


#plot ge associations
plot_grid(glvl_ge,
          gdiff_ge,
          nrow=2,
          labels = LETTERS[1:2])



# GE CHANGE AND GBM CHANGE -----------------------------------------------

#load the methylkit results for gbM differences
gene_mkit_list = list()
for (n in dat_names){
  print(n)
  mk_df = load_methylkit_results(n, 'gene_change_methylKit.Rdata')
  gene_mkit_list[[n]] = mk_df
}


#load the DESeq results files
ge_change_list = map(res_list, function(x) data.frame(x) %>% 
                       rownames_to_column('name') %>% 
                       dplyr::select(name, log2FoldChange, pvalue, padj))

#check the uploads' headers
head(ge_change_list[[n]])
head(gene_mkit_list[[n]])


#merge accross datasets
gbm_ge_change_list = list()
for (n in dat_names){
  print(n)
  gbm_ge_change_list[[n]] = merge_ge_meth_response(n, ge_change_list, gene_mkit_list) %>% 
    left_join(epi_state_list[[n]], by='name')
}


#function to plot responses against each other
plot_responses = function(n, diff_lvl_list, xcol, ycol, YLIM, XLIM){
  logo = logo_list[[n]]
  p_df = diff_lvl_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(p_df, xcol, ycol, '', '', ALPHA=0.1, ylim=YLIM, xlim=XLIM) +
    geom_hline(yintercept = 0, lty=2) +
    geom_smooth(method='lm', se=FALSE) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank())
  ggdraw(plt) + draw_image(logo,
                           x=logo_corner - logo_scale,
                           y=logo_corner-logo_scale, 
                           width = logo_scale,
                           height = logo_scale)
}


###### PLOT FOR ALL GENES

#build plots
gbm_response_plotlist = list()
for (n in dat_names){
  print(n)
  gbm_response_plotlist[[n]] = plot_responses(n, gbm_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}

#assemble multiplanel
ry = 1/30
rx = 1/15
gbm_response = assemble_row_panels(gbm_response_plotlist, 2, 'gbM difference', 'expression difference', ry, rx)
gbm_response


######### PLOT ALL GENES WITH VOLCANOS

volc_scatter_list = list()
for (n in dat_names){
  print(n)
  dat = gbm_ge_change_list[[n]]
  n_sig_ge = sum(dat$padj < 0.1, na.rm = TRUE)
  sig_ge_col = 'dodgerblue'
  ge_volc = plot_volcano_general(dat,
                                  xcol='log2FoldChange',
                                  ycol='pvalue.y',
                                  sigcol='padj',
                                  sigcut=0.1,
                                  ALPHA=0.5,
                                  title='',
                                  subtitle='',
                                  xlab='difference',
                                  ylab=bquote(log[10]*pvalue)) +
    scale_color_manual(values=c(sig_ge_col, 'black')) +
    theme(axis.title = element_blank(),
          legend.position = 'none',
          plot.subtitle = element_text(hjust = 0.5,
                                       color=sig_ge_col)) +
    labs(subtitle=paste('N=', n_sig_ge, sep='')) 
  n_sig_m = sum(dat$qvalue < 0.1, na.rm=TRUE)
  m_volc = plot_volcano_general(dat,
                                 xcol='meth.diff',
                                 ycol='pvalue.x',
                                 sigcol='qvalue',
                                 sigcut=0.1,
                                 ALPHA=0.5,
                                 title='',
                                 subtitle='',
                                 xlab='difference',
                                 ylab=bquote(log[10]*pvalue)) +
    theme(axis.title = element_blank(),
          legend.position = 'none',
          plot.subtitle = element_text(hjust = 0.5,
                                       color='red')) +
    labs(subtitle=paste('N=', n_sig_m, sep=''))
  
  volcs = plot_grid(m_volc, ge_volc, nrow=1)
  scatter = gbm_response_plotlist[[n]]
  plt = plot_grid(volcs, scatter, nrow=2)
  volc_scatter_list[[n]] = plt
}

plot_grid(plotlist = volc_scatter_list, nrow=2)


###### PLOT FOR METHYLATED GENES

#filter for methylated only
Mgbm_ge_change_list = map(gbm_ge_change_list, function(x) return(x %>% filter(gbm_class == 'methylated')))

#build plots
Mgbm_response_plotlist = list()
for (n in dat_names){
  print(n)
  Mgbm_response_plotlist[[n]] = plot_responses(n, Mgbm_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}

#assemble multiplanel
ry = 1/30
rx = 1/15
Mgbm_response = assemble_row_panels(Mgbm_response_plotlist, 2, 'GBM change', 'expression change', ry, rx)
Mgbm_response


###### PLOT FOR UNMETHYLATED GENES

#filter for methylated only
UMgbm_ge_change_list = map(gbm_ge_change_list, function(x) return(x %>% filter(gbm_class == 'unmethylated')))

#build plots
UMgbm_response_plotlist = list()
for (n in dat_names){
  print(n)
  UMgbm_response_plotlist[[n]] = plot_responses(n, UMgbm_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}

#assemble multiplanel
ry = 1/30
rx = 1/15
UMgbm_response = assemble_row_panels(UMgbm_response_plotlist, 2, 'GBM change', 'expression change', ry, rx)
UMgbm_response


# GE AND PROMOTER CHANGE --------------------------------------------------


#load the promoter methylation change results
promoter_mkit_list = list()
for (n in dat_names){
  print(n)
  mk_df = load_methylkit_results(n, 'promoter_change_methylKit.Rdata')
  promoter_mkit_list[[n]] = mk_df
}

#merge accross datasets
promoter_ge_change_list = list()
for (n in dat_names){
  print(n)
  promoter_ge_change_list[[n]] = merge_ge_meth_response(n, ge_change_list, promoter_mkit_list) %>% 
    left_join(epi_state_list[[n]], by='name')
}

##### PLOT ACCROSS ALL GENES

#build plots
promoter_response_plotlist = list()
for (n in dat_names){
  promoter_response_plotlist[[n]] = plot_responses(n, promoter_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}


#assemble multiplanel
ry = 1/30
rx = 1/15
promoter_response = assemble_row_panels(promoter_response_plotlist, 2, 'promoter methylation change', 'expression change', ry, rx)
promoter_response

######### PLOT ALL GENES WITH VOLCANOS

pvolc_scatter_list = list()
for (n in dat_names){
  print(n)
  dat = promoter_ge_change_list[[n]]
  n_sig_ge = sum(dat$padj < 0.1, na.rm = TRUE)
  sig_ge_col = 'dodgerblue'
  ge_volc = plot_volcano_general(dat,
                                 xcol='log2FoldChange',
                                 ycol='pvalue.y',
                                 sigcol='padj',
                                 sigcut=0.1,
                                 ALPHA=0.5,
                                 title='',
                                 subtitle='',
                                 xlab='difference',
                                 ylab=bquote(log[10]*pvalue)) +
    scale_color_manual(values=c(sig_ge_col, 'black')) +
    theme(axis.title = element_blank(),
          legend.position = 'none',
          plot.subtitle = element_text(hjust = 0.5,
                                       color=sig_ge_col)) +
    labs(subtitle=paste('N=', n_sig_ge, sep='')) 
  n_sig_m = sum(dat$qvalue < 0.1, na.rm=TRUE)
  m_volc = plot_volcano_general(dat,
                                xcol='meth.diff',
                                ycol='pvalue.x',
                                sigcol='qvalue',
                                sigcut=0.1,
                                ALPHA=0.5,
                                title='',
                                subtitle='',
                                xlab='difference',
                                ylab=bquote(log[10]*pvalue)) +
    theme(axis.title = element_blank(),
          legend.position = 'none',
          plot.subtitle = element_text(hjust = 0.5,
                                       color='red')) +
    labs(subtitle=paste('N=', n_sig_m, sep=''))
  
  volcs = plot_grid(m_volc, ge_volc, nrow=1)
  scatter = promoter_response_plotlist[[n]]
  plt = plot_grid(volcs, scatter, nrow=2)
  pvolc_scatter_list[[n]] = plt
}

plot_grid(plotlist = pvolc_scatter_list, nrow=2)




##### PLOT FOR METHYLATED ONLY

#filter for methylated only
Mpromoter_ge_change_list = map(promoter_ge_change_list, function(x) return(x %>% 
                                                                             filter(gbm_class=='methylated')))

#build plots
Mpromoter_response_plotlist = list()
for (n in dat_names){
  Mpromoter_response_plotlist[[n]] = plot_responses(n, Mpromoter_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}


#assemble multiplanel
ry = 1/30
rx = 1/15
Mpromoter_response = assemble_row_panels(Mpromoter_response_plotlist, 2, 'promoter methylation change', 'expression change', ry, rx)
Mpromoter_response


##### PLOT FOR Mbody Mpromoter only

#filter for methylated only
MgMp_promoter_ge_change_list = map(promoter_ge_change_list, function(x) return(x %>% 
                                                                             filter(epi_class=='Mg.Mp')))

#build plots
MgMp_promoter_response_plotlist = list()
for (n in dat_names){
  MgMp_promoter_response_plotlist[[n]] = plot_responses(n, MgMp_promoter_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}


#assemble multiplanel
ry = 1/30
rx = 1/15
MgMp_promoter_response = assemble_row_panels(MgMp_promoter_response_plotlist, 2, 'promoter methylation change', 'expression change', ry, rx)
MgMp_promoter_response


# PLOT PROMOTER LEVEL SEESAW ----------------------------------------------

##### FOR ALL GENES
pss_list = list()
for (n in dat_names){
  plt = plot_responses(n, promoter_ge_change_list, 'promoter_lvl', 'log2FoldChange', YLIM=FALSE, XLIM=FALSE)
  pss_list[[n]] = plt
}
allPss_plt = assemble_row_panels(pss_list, 2, 'promoter methylation level', 'expression change', ry, rx)
allPss_plt

##### FOR METHYLATED GENES
Mpss_list = list()
for (n in dat_names){
  plt = plot_responses(n, Mpromoter_ge_change_list, 'promoter_lvl', 'log2FoldChange', YLIM=FALSE, XLIM=FALSE)
  Mpss_list[[n]] = plt
}
Mpss_plt = assemble_row_panels(Mpss_list, 2, 'promoter methylation level', 'expression change', ry, rx)
Mpss_plt


##### FOR Mg_Mp genes
Mg_Mp_pss_list = list()
for (n in dat_names){
  plt = plot_responses(n, MgMp_promoter_ge_change_list, 'promoter_lvl', 'log2FoldChange', YLIM=FALSE, XLIM=FALSE)
  Mg_Mp_pss_list[[n]] = plt
}
Mg_Mp_pss_plt = assemble_row_panels(Mg_Mp_pss_list, 2, 'promoter methylation level', 'expression change', ry, rx)
Mg_Mp_pss_plt


# GE AND TSS CHANGE --------------------------------------------------


#load the promoter methylation change results
tss_mkit_list = list()
for (n in dat_names){
  print(n)
  mk_df = load_methylkit_results(n, 'tss_change_methylKit.Rdata')
  tss_mkit_list[[n]] = mk_df
}

#merge accross datasets
tss_ge_change_list = list()
for (n in dat_names){
  print(n)
  tss_ge_change_list[[n]] = merge_ge_meth_response(n, ge_change_list, tss_mkit_list)
}


#build plots
tss_response_plotlist = list()
for (n in dat_names){
  promoter_response_plotlist[[n]] = plot_responses(n, tss_ge_change_list, 'meth.diff', 'log2FoldChange', YLIM=c(-3,3), XLIM=c(-10,10))
}


#assemble multiplanel
ry = 1/30
rx = 1/15
tss_response_plotlist = assemble_row_panels(promoter_response_plotlist, 2, 'TSS methylation change', 'expression change', ry, rx)


# LOOK AT SEESAW ----------------------------------------------------------


diff_lvl_list = list()
for (n in dat_names){
  diff_lvl_list[[n]] = merge_ge_meth_response(n, gene_lvl_list, gene_mkit_list)
}

#plot seesaw
plot_seesaw = function(n, diff_lvl_list, xcol, ycol, YLIM, logo.x=0.25, logo.y=0.65){
  logo = logo_list[[n]]
  p_df = diff_lvl_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(p_df, xcol, ycol, '', '', ALPHA=0.1, ylim=YLIM) +
    geom_hline(yintercept = 0, lty=2) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank())
  ggdraw(plt) + draw_image(logo,
                           x=logo.x,
                           y=logo.y,
                           width = logo_scale,
                           height = logo_scale)
}

seesaw_plotlist = list()
for (n in dat_names){
  seesaw_plotlist[[n]] = plot_seesaw(n, diff_lvl_list, 'lfrac_meth', 'meth.diff', c(-20, 20))
}

#assemble multiplanel
ry = 1/30
rx = 1/15
m_seesaw = assemble_row_panels(seesaw_plotlist, 2, 'gbM level', 'gbM difference', ry, rx)
m_seesaw


# CHANGE DISTRIBUTIONS FOR METHYLATION CLASSES ----------------------------
#idea here is to look at how each methylation class changes
#for a true seesaw, each class should be shifted oppositely from zero

dat_names
select_dats = c('silkworm', 'termite', 'maternal')


plot_meth_class_density = function(n, diff_lvl_list, difference_col, add_logo=TRUE, xlab=''){
  ss_df = diff_lvl_list[[n]]
  logo = logo_list[[n]]
  ss_df %>% 
    ggplot(aes(x=lfrac_meth)) +
    geom_histogram() + 
    geom_vline(xintercept = -6) +
    scale_x_continuous(labels = log2_to_percent)
  
  
  plt=ss_df %>% 
    mutate(meth_class = if_else(lfrac_meth > -6,
                                'methylated',
                                'unmethylated'),
           meth_class = factor(meth_class, levels=c('unmethylated', 'methylated'))) %>% 
    ggplot(aes_string(x=difference_col)) +
    geom_density() + 
    geom_vline(xintercept = 0, lty=2) +
    labs(x=xlab) +
    facet_wrap(~meth_class, scales='free')
  if (add_logo){
    plt2 = ggdraw(plt) + draw_image(logo,
                             x=logo_corner - logo_scale + 0.06,
                             y=logo_corner-logo_scale-0.1, 
                             width = logo_scale,
                             height = logo_scale)
  } else {
    plt2=plt
  }
    return(plt2)
}
  
#build plots accross select datasets
seesaw_dens_list = list()
for (n in select_dats){
  seesaw_dens_list[[n]] = plot_meth_class_density(n, diff_lvl_list, 'meth.diff', xlab='gbM difference')
}

#plot
mss_dens_plt = plot_grid(plotlist = seesaw_dens_list,
                         nrow=length(select_dats))


# PLOT GE SEESAW ----------------------------------------------------------
#repeat above for gene expression



#merge with gbm level
ge_diff_lvl_list = list()
for (n in dat_names){
  print(n)
  ge_diff_lvl_list[[n]] = merge_ge_meth_response(n, gene_lvl_list, ge_change_list)
}

#build plot
ge_seesaw_plotlist = list()
for (n in dat_names){
  ge_seesaw_plotlist[[n]] = plot_seesaw(n,
                                        ge_diff_lvl_list,
                                        'lfrac_meth',
                                        'log2FoldChange',
                                        c(-3, 3))
}

#assemble multiplanel
ry = 1/30
rx = 1/15
ge_seesaw = assemble_row_panels(ge_seesaw_plotlist, 2, 'gbM level', 'expression difference', ry, rx)
ge_seesaw

# GE CHANGE DISTRIBUTIONS FOR METHYLATION CLASSES ----------------------------


dat_names
select_dats = c('bumblebee', 'subcaste')


#build plots accross select datasets
ge_seesaw_dens_list = list()
for (n in select_dats){
  print(n)
  ge_seesaw_dens_list[[n]] = plot_meth_class_density(n, ge_diff_lvl_list, 'log2FoldChange', add_logo=FALSE) +
    labs(x='expression difference') +
    lims(x=c(-2,2)) +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

#plot
gss_dens_pans = plot_grid(plotlist = ge_seesaw_dens_list,
                         nrow=length(select_dats))
xlab = ggdraw() + draw_label(' ')
gss_dens_plt = plot_grid(gss_dens_pans, xlab, nrow=2, rel_heights = c(20,1))


# assemble ----------------------------------------------------------------

plot_sub_seesaw = function(n, diff_lvl_list, xcol, ycol, YLIM, logo.x=0.25, logo.y=0.65){
  logo = logo_list[[n]]
  p_df = diff_lvl_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(p_df, xcol, ycol, '', '', ALPHA=0.1, ylim=YLIM) +
    geom_hline(yintercept = 0, lty=2) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(labels = log2_to_percent)
  ggdraw(plt) + draw_image(logo,
                           x=logo.x,
                           y=logo.y,
                           width = logo_scale,
                           height = logo_scale)
}

sub_seesaws = list()
for (n in select_dats){
  sss = plot_sub_seesaw(n,
                    ge_diff_lvl_list,
                    'lfrac_meth',
                    'log2FoldChange',
                    c(-1, 1),
                    logo.x = 0.7,
                    logo.y = 0.8)
  sub_seesaws[[n]] = sss
}

ss_pans = plot_grid(plotlist = sub_seesaws, nrow=2)
plot_grid(ss_pans, gss_dens_plt, nrow=1, rel_widths = c(.75, 1))

#COMPUTE THE VARIANCE IN DIFFERENTIAL EXPRESSION EXPLAINED BY THE GROUPS

dat = ge_diff_lvl_list[['bumblebee']] %>% 
  mutate(meth_class = if_else(lfrac_meth > -6,
                              'methylated',
                              'unmethylated'),
         meth_class = factor(meth_class, levels=c('unmethylated', 'methylated'))) %>% 
  as_tibble()

lm1 = lm(log2FoldChange ~ meth_class, data=dat)
summary(lm1)
aov(lm1)


# GE CHANGE BY METHYLATION CLASS ------------------------------------------

n='silkworm'
df = ge_diff_lvl_list[[n]]

df %>% 
  mutate(meth_class = if_else(lfrac_meth > -6,
                              'methylated',
                              'unmethylated'),
         meth_class = factor(meth_class, levels=c('unmethylated', 'methylated'))) %>% 
  ggplot(aes(x=lfrac_meth, y=log2FoldChange)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_point(alpha=0.1) +
  geom_smooth(aes(color=meth_class), method='lm')





# PLOT GE VS INDIVIDUAL EXON ----------------------------------------------

#load the exon methylation change results
exon_mkit_list = list()
for (n in dat_names){
  print(n)
  mk_df = load_methylkit_results(n, 'exon_change_methylKit.Rdata')
  bpd = bioprojectDirs[n]
  cds_gene_file = paste('bioprojects/', bpd, '/cdsID_geneID.tsv', sep='')
  cds_gene = read_tsv(cds_gene_file) %>% 
    dplyr::rename(name=cds)
  exon_df = mk_df %>% 
    left_join(cds_gene, by = 'name') %>% 
    group_by(gene) %>% 
    summarize(meth.diff = max(meth.diff)) %>% 
    ungroup() %>% 
    dplyr::rename(name=gene)
  exon_mkit_list[[n]] = exon_df
}


#merge with ge
exon_ge_list = list()
for (n in dat_names){
  print(n)
  gedat = ge_change_list[[n]]
  edat = exon_mkit_list[[n]]
  exon_ge_list[[n]] = full_join(edat, gedat, by = 'name') %>% 
    mutate(absl2 = abs(log2FoldChange),
           abs.meth.diff = abs(meth.diff))
}


#build plots
e_plt_list = list()
eabs_plt_list = list()
for (n in dat_names){
  e_plt_list[[n]] = plot_responses(n, exon_ge_list, 'meth.diff', 'log2FoldChange', YLIM=FALSE)
  eabs_plt_list[[n]] = plot_responses(n, exon_ge_list, 'abs.meth.diff', 'absl2', YLIM=FALSE)
}

ae_plts = assemble_row_panels(e_plt_list, 2, 'max exon methylation change', 'expression change', ry, rx)
aeabs_plts = assemble_row_panels(eabs_plt_list, 2, 'max exon methylation change', 'abs expression change', ry, rx)

# CHANGES IN WINDOWS ------------------------------------------------------

#function to read in window meth data and merge with ge
read_windows = function(n, window_file, dist_cut){
  bpd = bioprojectDirs[n]
  window_path = paste('bioprojects/', bpd, '/meth_response/significant/', window_file, sep='')
  ge_diffs = ge_change_list[[n]]
  w_df = read_tsv(window_path)
  if (nrow(w_df) > 0){
    res = w_df %>% 
      mutate(split1 = sapply(gene.description, function(x) strsplit(x, 'ID=')[[1]][2]),
             name = sapply(split1, function(x) strsplit(x, ';')[[1]][1])) %>% 
      dplyr::select(-gene.chr, -gene.description, -gene.strand, -split1)
      
  }
  else {
    res = 'no significant'
  }
}


window_file = '1KbWindows_change_significant_closeGeneDistances.tsv'
min_dist = 2000
require_in_gene = TRUE


#upload the 500bp windows
window_list0 = list()
for (n in dat_names){
  print(n)
  res = read_windows(n, '500bpWindows_change_significant_closeGeneDistances.tsv', 2e3)
  if (class(res)=="character"){
    next
  } else{
    window_list0[[n]] = res
  }
}
names(window_list0)


filter_windows = function(x){
  if (require_in_gene){
    x= x %>% 
      filter(gene.distance < min_dist,
             gene.enclosesWindow)
  } else {
    x = x %>% 
      filter(gene.distance < min_dist)
  }
  return(x)
}

window_list = map(window_list0, function(x) filter_windows(x))


make_ge_calls = function(n){
  print(n)
  n=dat_names[1]
  wdat = window_list[[n]]
  sig_meth = wdat$name
  gedat = ge_change_list[[n]] %>% 
    mutate(isSig = name %in% sig_meth,
           absl2 = abs(log2FoldChange)) %>% 
    as_tibble()
  return(gedat)
}

ge_windows = map(dat_names, make_ge_calls)
names(ge_windows) = dat_names


#function to plot l2
plot_l2 = function(n){
  gedat = ge_windows[[n]]
  logo = logo_list[[n]]
  plt = gedat %>% 
    ggplot(aes(x=isSig, y=log2FoldChange)) +
    geom_boxplot()
  ggdraw(plt) + draw_image(logo,
                           x=logo_corner - logo_scale,
                           y=logo_corner-logo_scale, 
                           width = logo_scale,
                           height = logo_scale)
}

l2_plts = map(dat_names, plot_l2)
l2_ass = plot_grid(plotlist = l2_plts, nrow=2)


plot_absl2 = function(n){
  gedat = ge_windows[[n]]
  logo = logo_list[[n]]
  plt = gedat %>% 
    ggplot(aes(x=isSig, y=absl2)) +
    geom_boxplot() +
    lims(y=c(0,1))
  ggdraw(plt) + draw_image(logo,
                           x=logo_corner - logo_scale,
                           y=logo_corner-logo_scale, 
                           width = logo_scale,
                           height = logo_scale)
}

absl2_plts = map(dat_names, plot_absl2)
absl2_ass = plot_grid(plotlist = absl2_plts, nrow=2)

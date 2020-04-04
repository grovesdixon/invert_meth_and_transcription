#plot_meth_level_figs.R

rm(list=ls())
source('invert_meth_ge_functions.R')
source('metadata/species_logos/logo_source.R')


#GLOBAL VARS
bioprojectDirs


# LOAD THE GBM LEVEL DATA ---------------------------------------------------

#set up the path list
gene_path_list = paste('bioprojects/', bioprojectDirs, '/meth_level/gene_basicStatsBed.tsv', sep='')
exon_path_list = paste('bioprojects/', bioprojectDirs, '/meth_level/exon_basicStatsBed.tsv', sep='')
names(gene_path_list) = dat_names
names(exon_path_list) = dat_names

#funciton to read in rld data
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
    mutate(lfrac_meth = log(fracMeth, 2)) %>% 
    dplyr::select(-fracMeth)
}

#get the rld data for each project
gene_lvl_list = map(gene_path_list, load_lvl)
exon_lvl_list0 = map(exon_path_list, load_lvl)


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
  ghist = plot_bimodal(gdf, logo, logo_scale)
  ehist = plot_bimodal(edf, logo, logo_scale)
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
    scale_x_continuous(labels = log2_to_percent) +
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
for (n in dat_names){
  glvls_plt_list[[n]] = plot_lvl_lvl(n, gene_lvl_list)
  elvls_plt_list[[n]] = plot_lvl_lvl(n, exon_lvl_list)
}


#assmeble multiplanels
ry = 1/30
rx = 1/15
glvl_ge = assemble_row_panels(glvls_plt_list, 2, '% methylation', 'expression level', ry, rx)
elvl_ge = assemble_row_panels(elvls_plt_list, 2, '% methylation', 'expression level', ry, rx)


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
    scale_x_continuous(labels = log2_to_percent) +
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
glvl_ge = assemble_row_panels(g_change_plt_list, 2, '% methylation', 'expression difference', ry, rx)
elvl_ge = assemble_row_panels(e_change_plt_list, 2, '% methylation', 'expression difference', ry, rx)



# LOOK AT SEESAW ----------------------------------------------------------

#function to load methylation differences from methylkit
load_methylkit_results = function(n, file_name){
  bpd = bioprojectDirs[n]
  in_path = paste('bioprojects/', bpd, '/meth_response/', file_name, sep='')
  ll = load(in_path)
  df = data.frame(myDiff) %>% 
    left_join(wbounds, by = c('chr', 'start', 'end')) %>% 
    dplyr::select(name, pvalue, qvalue, meth.diff) %>% 
    as_tibble()
}

#load the results
gene_mkit_list = list()
for (n in dat_names){
  mk_df = load_methylkit_results(n, 'gene_change_methylKit.Rdata')
  gene_mkit_list[[n]] = mk_df
}


#merge with level data
merge_diff_lvl = function(n, change_list){
  lvl_df = gene_lvl_list[[n]]
  diff_df = change_list[[n]]
  m_df = inner_join(lvl_df, diff_df, by = 'name')
}

diff_lvl_list = list()
for (n in dat_names){
  diff_lvl_list[[n]] = merge_diff_lvl(n, gene_mkit_list)
}

#plot seesaw
plot_seesaw = function(n, diff_lvl_list, xcol, ycol, YLIM){
  logo = logo_list[[n]]
  p_df = diff_lvl_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(p_df, xcol, ycol, '', '', ALPHA=0.1, ylim=YLIM) +
    geom_hline(yintercept = 0, lty=2) +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_continuous(labels = log2_to_percent) +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(axis.title = element_blank())
  ggdraw(plt) + draw_image(logo,
                           x=0.2,
                           y=0.65,
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
m_seesaw = assemble_row_panels(seesaw_plotlist, 2, '% methylation', 'methylation change', ry, rx)



# CHANGE DISTRIBUTIONS FOR METHYLATION CLASSES ----------------------------
#idea here is to look at how each methylation class changes
#for a true seesaw, each class should be shifted oppositely from zero

dat_names
select_dats = c('silkworm', 'termite', 'maternal')


plot_meth_class_density = function(n, diff_lvl_list, difference_col){
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
    facet_wrap(~meth_class, scales='free')
    plt2 = ggdraw(plt) + draw_image(logo,
                             x=logo_corner - logo_scale + 0.06,
                             y=logo_corner-logo_scale-0.1, 
                             width = logo_scale,
                             height = logo_scale)
    return(plt2)
}
  
#build plots accross select datasets
seesaw_dens_list = list()
for (n in select_dats){
  seesaw_dens_list[[n]] = plot_meth_class_density(n, diff_lvl_list, 'meth.diff')
}

#plot
mss_dens_plt = plot_grid(plotlist = seesaw_dens_list,
                         nrow=length(select_dats))


# PLOT GE SEESAW ----------------------------------------------------------
#repeat above for gene expression

#reformat ge results
ge_change_list = map(res_list, function(x) data.frame(x) %>% 
                       rownames_to_column('name') %>% 
                       dplyr::select(name, log2FoldChange, pvalue, padj))

#merge with gbm level
ge_diff_lvl_list = list()
for (n in dat_names){
  print(n)
  ge_diff_lvl_list[[n]] = merge_diff_lvl(n, ge_change_list)
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
ge_seesaw = assemble_row_panels(ge_seesaw_plotlist, 2, '% methylation', 'expression change', ry, rx)


# GE CHANGE DISTRIBUTIONS FOR METHYLATION CLASSES ----------------------------


dat_names
select_dats = c('bumblebee', 'subcaste')


#build plots accross select datasets
ge_seesaw_dens_list = list()
for (n in select_dats){
  print(n)
  ge_seesaw_dens_list[[n]] = plot_meth_class_density(n, ge_diff_lvl_list, 'log2FoldChange')
}

#plot
gss_dens_plt = plot_grid(plotlist = ge_seesaw_dens_list,
                         nrow=length(select_dats))



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


# GE AND GBM CHANGE -----------------------------------------------

#function to merge meth response and ge response dfs
head(ge_change_list[[n]])
head(gene_mkit_list[[n]])

#function to merge
merge_ge_meth_response = function(n, ge_list, meth_list){
  mdf = meth_list[[n]]
  gdf = ge_list[[n]]
  inner_join(mdf, gdf, by='name')
}

#merge accross datasets
gbm_ge_change_list = list()
for (n in dat_names){
  print(n)
  gbm_ge_change_list[[n]] = merge_ge_meth_response(n, ge_change_list, gene_mkit_list)
}


#function to plot responses against each other
plot_responses = function(n, diff_lvl_list, xcol, ycol, YLIM, XLIM){
  logo = logo_list[[n]]
  p_df = diff_lvl_list[[n]]
  plt = plot_scatter_pearsonCor_annotated(p_df, xcol, ycol, '', '', ALPHA=0.1, ylim=YLIM) +
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


#build plots
gbm_response_plotlist = list()
for (n in dat_names){
  print(n)
  gbm_response_plotlist[[n]] = plot_responses(n, gbm_ge_change_list, 'meth.diff', 'log2FoldChange', FALSE)
}


#assemble multiplanel
ry = 1/30
rx = 1/15
gbm_response = assemble_row_panels(gbm_response_plotlist, 2, 'GBM change', 'expression change', ry, rx)



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
  promoter_ge_change_list[[n]] = merge_ge_meth_response(n, ge_change_list, promoter_mkit_list)
}


#build plots
promoter_response_plotlist = list()
for (n in dat_names){
  promoter_response_plotlist[[n]] = plot_responses(n, promoter_ge_change_list, 'meth.diff', 'log2FoldChange', FALSE)
}


#assemble multiplanel
ry = 1/30
rx = 1/15
promoter_response = assemble_row_panels(promoter_response_plotlist, 2, 'promoter methylation change', 'expression change', ry, rx)


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
  promoter_response_plotlist[[n]] = plot_responses(n, tss_ge_change_list, 'meth.diff', 'log2FoldChange', FALSE)
}


#assemble multiplanel
ry = 1/30
rx = 1/15
tss_response_plotlist = assemble_row_panels(promoter_response_plotlist, 2, 'TSS methylation change', 'expression change', ry, rx)



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

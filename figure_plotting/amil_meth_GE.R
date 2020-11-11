#amil_meth_ge.R
rm(list=ls())

source('invert_meth_ge_functions.R')

# LOAD GBM DATA -----------------------------------------------------------

ll=load('bioprojects/amillepora_PRJNA601565/meth_level/gbmLvl.Rdata')
ll
gbm.dat = gbm.dat %>% 
  dplyr::select(name, l.fracMeth, mrB, mbd.score) %>% 
  as_tibble()

ll=load('bioprojects/amillepora_PRJNA601565/meth_level/exonGbmLvl.Rdata')
ll

eg.dat=eg.dat %>% 
  dplyr::rename(l.fracMeth=lfracMeth) %>% 
  mutate(l.fracMeth = if_else(l.fracMeth < SET_ZERO,
                              SET_ZERO,
                              l.fracMeth))


# LOAD EXPRESSION DATA --------------------------------------------------

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/rld.Rdata')
ll

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/genotype_results.Rdata')
ll
g.res=res %>% 
  data.frame() %>% 
  mutate(name=rownames(res)) %>% 
  dplyr::select(name, log2FoldChange, pvalue, padj) %>% 
  as_tibble()

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/tissue_results.Rdata')
ll
t.res=res %>% 
  data.frame() %>% 
  mutate(name=rownames(res)) %>% 
  dplyr::select(name, log2FoldChange, pvalue, padj) %>% 
  as_tibble()



# PLOT HISTOGRAMS ---------------------------------------------------------

#wgbs
Npm = sum(!is.na(eg.dat$l.fracMeth))
PM.BREAKS = log(c(2^SET_ZERO, .02, .15,1), 2)
PM.XLIM = c(PM.BREAKS[1], PM.BREAKS[length(PM.BREAKS)])
pmHist = eg.dat %>% 
  ggplot(aes(x=l.fracMeth)) +
  geom_histogram() +
  scale_x_continuous(breaks = PM.BREAKS,
                     labels = log2_to_percent) +
  labs(title='WGBS', subtitle=paste(Npm, 'genes'), x='% methylation')

#mbdseq
Nmbd = sum(!is.na(gbm.dat$mbd.score))
MBD.BREAKS = c(-2, 0, 2, 4)
MBD.XLIM = c(MBD.BREAKS[1], MBD.BREAKS[length(MBD.BREAKS)])
mbdHist = gbm.dat %>% 
  ggplot(aes(x=mbd.score)) +
  geom_histogram() +
  scale_x_continuous(breaks = MBD.BREAKS,
                     limits = c(MBD.BREAKS[1], MBD.BREAKS[length(MBD.BREAKS)])) +
  labs(title='MBD-seq', subtitle=paste(Nmbd, 'genes'),x='MBD-score')

#mdRAD
Nmr = sum(!is.na(gbm.dat$mrB))
MR.BREAKS = c(-8, 0, 8)
MR.XLIM = c(MR.BREAKS[1], MR.BREAKS[length(MR.BREAKS)])
mrHist = gbm.dat %>% 
  ggplot(aes(x=mrB)) +
  geom_histogram() +
  scale_x_continuous(breaks = MR.BREAKS,
                     limits = c(MR.BREAKS[1], MR.BREAKS[length(MR.BREAKS)])) +
  labs(title='mdRAD', subtitle=paste(Nmr, 'genes'), x=bquote(log[2]~FPKM))

histList = list(pmHist, mbdHist, mrHist)

# COMPARE WITH EXPRESION LEVEL --------------------------------------------

#for genes
rld.df = data.frame(assay(rld))
genes = rownames(rld.df)
ge.lvl = apply(rld.df, 1, function(x) mean(x, na.rm=TRUE)) %>% 
  data.frame() %>% 
  mutate(name=genes) %>% 
  as_tibble() %>% 
  set_names('mnGe', 'name')

#merge with expression level
gbm.ge = gbm.dat %>% 
  left_join(ge.lvl, by = 'name') 

#WGBS for genes
notUsing = plot_scatter_pearsonCor_annotated(gbm.ge, 'l.fracMeth', 'mnGe', 'GBM', 'expression level', ALPHA=0.1) +
  geom_smooth(method='lm', se=FALSE)

#mbd-seq genes
mbdLvl = plot_scatter_pearsonCor_annotated(gbm.ge, 'mbd.score', 'mnGe', 'MBD-score GBM', 'expression level',
                                           ALPHA=0.1, xlim=MBD.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

#methylRAD for genes
radLvl = plot_scatter_pearsonCor_annotated(gbm.ge, 'mrB', 'mnGe', 'mdRAD GBM', 'expression level',
                                           ALPHA=0.1, xlim=MR.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5))


#for exons
eg.dat
e.ge = eg.dat %>% 
  left_join(ge.lvl, by = 'name')

wgbsLvl = plot_scatter_pearsonCor_annotated(e.ge, 'l.fracMeth', 'mnGe', 'GBM', 'expression level',
                                            ALPHA=0.1, xlim=PM.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))


#mbd-exons
eMbdLvl=plot_scatter_pearsonCor_annotated(e.ge, 'mbd.score', 'mnGe', 'MBD-seq GBM', 'expression level',
                                          ALPHA=0.1, xlim=MBD.XLIM) +
  geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) 

#mbd-exons
eMdRADLvl=plot_scatter_pearsonCor_annotated(e.ge, 'mrB', 'mnGe', 'mdRAD GBM', 'expression level',
                                            ALPHA=0.1, xlim=MR.XLIM) +
  geom_smooth(method='lm', se=FALSE) 


#PLOT BEST TOGETHER
lvlList = list(wgbsLvl + scale_x_continuous(breaks = PM.BREAKS,
                                            labels = log2_to_percent),
               mbdLvl + scale_x_continuous(breaks = MBD.BREAKS,
                                           limits = c(MBD.BREAKS[1], MBD.BREAKS[length(MBD.BREAKS)])),
               radLvl + scale_x_continuous(breaks = MR.BREAKS,
                                           limits = c(MR.BREAKS[1], MR.BREAKS[length(MR.BREAKS)])))

# COMPARE WITH DIFFERENTIAL EXPRESSION ------------------------------------

YLIM = c(0,1)

gTissue = t.res %>% 
  left_join(gbm.dat, by='name') %>% 
  mutate(absl = abs(log2FoldChange))

eTissue = t.res %>% 
  left_join(eg.dat, by='name') %>% 
  mutate(absl = abs(log2FoldChange))
  



#plot for wgbs
wgbs.t = plot_scatter_pearsonCor_annotated(gTissue, 'l.fracMeth', 'absl', 'WGBS', 'tissue difference',
                                           ALPHA=0.1, ylim=YLIM, xlim=PM.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_x_continuous(breaks = PM.BREAKS,
                     labels = log2_to_percent)

#plot for mbd-seq
mbd.t = plot_scatter_pearsonCor_annotated(gTissue, 'mbd.score', 'absl', 'MBD-seq GBM', 'tissue difference',
                                          ALPHA=0.1, ylim=YLIM, xlim = MBD.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

#plot for methylRAD
mr.t = plot_scatter_pearsonCor_annotated(gTissue, 'mrB', 'absl', 'MethylRAD GBM', 'tissue difference',
                                         ALPHA=0.1, ylim=YLIM, xlim=MR.XLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5))

#plot for wgbs exons
wgbs.t=plot_scatter_pearsonCor_annotated(eTissue, 'l.fracMeth', 'absl', 'WGBS GBM', 'tissue difference', ALPHA=0.1, ylim=YLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_x_continuous(labels = log2_to_percent)

#plot for mbeseq exons
eMbdTip=plot_scatter_pearsonCor_annotated(eTissue, 'mbd.score', 'absl', 'MBD-score GBM', 'tissue difference', ALPHA=0.1, ylim=YLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2)

#plot for mdRAD
eMdTip = plot_scatter_pearsonCor_annotated(eTissue, 'mrB', 'absl', 'mdRAD GBM', 'tissue difference', ALPHA=0.1, ylim=YLIM) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2)


#plot together 
tipList = list(wgbs.t + scale_x_continuous(breaks = PM.BREAKS,
                                           labels = log2_to_percent),
               mbd.t + scale_x_continuous(breaks = MBD.BREAKS,
                                            limits = c(min(MBD.BREAKS), max(MBD.BREAKS))),
               mr.t + scale_x_continuous(breaks = MR.BREAKS,
                                         limits = c(min(MR.BREAKS), max(MR.BREAKS))))
# plot_grid(plotlist = tipList, nrow=1)



# ASSEMBLE TOGETHER -------------------------------------------------------


histList2 = lapply(histList, function(x) return(x+theme(axis.title = element_blank(),
                                                        plot.title = element_text(hjust=0.5),
                                                        plot.subtitle = element_text(hjust=0.5))))
lvlList2 = lapply(lvlList, function(x) return(x+theme(axis.title = element_blank(),
                                                      plot.title = element_blank()) +
                                                scale_y_continuous(labels=scaleFUN2)))
tipList2 = lapply(tipList, function(x) return(x + 
                                                theme(axis.title = element_blank(),
                                                      plot.title = element_blank()) +
                                                scale_y_continuous(limits = c(0,1),
                                                                   breaks = c(0, .25, .5, .75, 1),
                                                                   labels = c('0.00', '0.25', '0.50', '0.75', ''))))

hists = plot_grid(plotlist = histList2, nrow=1, labels = LETTERS[1:3])
lvls = plot_grid(plotlist = lvlList2, nrow=1, labels = LETTERS[4:6])
tips = plot_grid(plotlist = tipList2, nrow=1, labels = LETTERS[7:9])

hy =  ggdraw() + draw_label('count', angle=90)
ly =  ggdraw() + draw_label('mRNA level', angle=90)
ty =  ggdraw() + draw_label('mRNA difference', angle=90)
ywidth=1/30
yhists = plot_grid(hy, hists, nrow = 1, rel_widths = c(ywidth,1))
ylvls = plot_grid(ly, lvls, nrow = 1, rel_widths = c(ywidth,1))
ytips = plot_grid(ty, tips, nrow = 1, rel_widths = c(ywidth,1))

plts = plot_grid(yhists,
          ylvls,
          ytips,
          nrow=3)
xlab = ggdraw() + draw_label('gbM level')
plot_grid(plts, xlab, nrow=2, rel_heights = c(1,1/20))



# DIFFERENTIAL METH AND GE FOR GBM-----------------------------------

#LOAD GBM RESULTS FROM BENCHMARKING STUDY
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/tissue_gbm_response.Rdata')
ll  
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/genotype_gbm_response.Rdata')
ll

#merge with gene expression results
merge_meth_and_ge = function(gbm.df, res){
  colnames(res) = c('name', 'ge.log2FoldChange', 'ge.pvalue', 'ge.padj')
  gbm.df %>% 
    left_join(res, by = 'name') %>% 
    as_tibble()
    
}


#pick the contrast to look at
ge.gbm = merge_meth_and_ge(g.gbm.dat, g.res) #genotype (big differences, only mdRAD correlates)
ge.gbm = merge_meth_and_ge(t.gbm.dat, t.res) #tissue

##### BUILD VOLCANOS
#WGBS
n_sig_wgbs=sum(ge.gbm$qvalue < 0.1, na.rm=TRUE)
wgbs_volc = plot_volcano_general(ge.gbm,
                                 xcol='meth.diff',
                                 ycol='pvalue',
                                 sigcol='qvalue',
                                 xlab='difference',
                                 ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('red', 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_wgbs, sep=''),
       title = 'WGBS')

#MBD
n_sig_mbd=sum(ge.gbm$mbd.padj < 0.1, na.rm=TRUE)
mbd_volc = plot_volcano_general(ge.gbm,
                                 xcol='mbd.log2FoldChange',
                                 ycol='mbd.pvalue',
                                 sigcol='mbd.padj',
                                 xlab='difference',
                                 ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('red', 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_mbd, sep=''),
       title = 'MBD')

#mdRAD
n_sig_mdr=sum(ge.gbm$mr.padj < 0.1, na.rm=TRUE)
mr_volc = plot_volcano_general(ge.gbm,
                                xcol='mr.log2FoldChange',
                                ycol='mr.pvalue',
                                sigcol='mr.padj',
                                xlab='difference',
                                ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('red', 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_mdr, sep=''),
       title = 'mdRAD')

#GE
sig_ge_col = 'dodgerblue'
n_sig_ge=sum(ge.gbm$ge.padj < 0.1, na.rm=TRUE)
ge_volc = plot_volcano_general(ge.gbm,
                                xcol='ge.log2FoldChange',
                                ycol='ge.pvalue',
                                sigcol='ge.padj',
                                xlab='difference',
                                ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c(sig_ge_col, 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color=sig_ge_col)) +
  labs(subtitle=paste('N=', n_sig_ge, sep=''),
       title = 'GE')

#PLOT TOGETHER
volc_list = list(wgbs_volc, mbd_volc, mr_volc, ge_volc)
volc_plts = plot_grid(plotlist = volc_list, nrow=1)


##### BUILD SCATTERPLOTS

wgbs.scat = plot_scatter_pearsonCor_annotated(ge.gbm, 'meth.diff', 'ge.log2FoldChange', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)

mbd.scat = plot_scatter_pearsonCor_annotated(ge.gbm, 'mbd.log2FoldChange', 'ge.log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)


mr.scat = plot_scatter_pearsonCor_annotated(ge.gbm, 'mr.log2FoldChange', 'ge.log2FoldChange', 'mdRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)


#plot volcanos and scatterplots
scatter_list = list(wgbs.scat, mbd.scat, mr.scat)
mod_scatters = lapply(scatter_list, function(x) return(x + theme(axis.title = element_blank())))
scatter_plts = plot_grid(plotlist = mod_scatters, nrow=1)
vylab = ggdraw() + draw_label(bquote('-log'[10]~'p-value'), angle=90)
vplts = plot_grid(vylab, volc_plts, rel_widths = c(1,10))
sylab = ggdraw() + draw_label('expression difference', angle=90)
splts = plot_grid(sylab, scatter_plts, rel_widths = c(1,12))
xlab = ggdraw() + draw_label('gbM difference')
# plot_grid(vplts, splts, xlab, nrow = 3, rel_heights = c(12,12,1))


#or plot with just scatterplots
plot_grid(splts, xlab, nrow=2, rel_heights = c(8,1))

# FOR PROMOTERS -----------------------------------------------------------

#LOAD PROMOTER DIFFERENCES FROM BENCHMARKING STUDY
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/tissue_promoter_response.Rdata')
ll  
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/genotype_promoter_response.Rdata')
ll  

#pick contrast to plot
ge.prom = merge_meth_and_ge(g.dat, g.res) #genotype
ge.prom = merge_meth_and_ge(t.dat, t.res) #tissue

##### BUILD VOLCANOS
#WGBS
n_sig_wgbs=sum(ge.prom$qvalue < 0.1, na.rm=TRUE)
wgbs_volc = plot_volcano_general(ge.prom,
                                 xcol='meth.diff',
                                 ycol='pvalue',
                                 sigcol='qvalue',
                                 xlab='difference',
                                 ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('red', 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_wgbs, sep=''),
       title = 'WGBS')

#MBD
n_sig_mbd=sum(ge.prom$mbd.padj < 0.1, na.rm=TRUE)
mbd_volc = plot_volcano_general(ge.prom,
                                xcol='mbd.log2FoldChange',
                                ycol='mbd.pvalue',
                                sigcol='mbd.padj',
                                xlab='difference',
                                ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('red', 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_mbd, sep=''),
       title = 'MBD')

#mdRAD
n_sig_mdr=sum(ge.prom$mr.padj < 0.1, na.rm=TRUE)
mr_volc = plot_volcano_general(ge.prom,
                               xcol='mr.log2FoldChange',
                               ycol='mr.pvalue',
                               sigcol='mr.padj',
                               xlab='difference',
                               ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c('black', 'red')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color='red')) +
  labs(subtitle=paste('N=', n_sig_mdr, sep=''),
       title = 'mdRAD')

#GE
sig_ge_col = 'dodgerblue'
n_sig_ge=sum(ge.prom$ge.padj < 0.1, na.rm=TRUE)
ge_volc = plot_volcano_general(ge.prom,
                               xcol='ge.log2FoldChange',
                               ycol='ge.pvalue',
                               sigcol='ge.padj',
                               xlab='difference',
                               ylab=bquote(log[10]*pvalue)) +
  scale_color_manual(values=c(sig_ge_col, 'black')) +
  theme(axis.title = element_blank(),
        legend.position = 'none',
        plot.subtitle = element_text(hjust = 0.5,
                                     color=sig_ge_col)) +
  labs(subtitle=paste('N=', n_sig_ge, sep=''),
       title = 'GE')

#PLOT TOGETHER
volc_list = list(wgbs_volc, mbd_volc, mr_volc, ge_volc)
volc_plts = plot_grid(plotlist = volc_list, nrow=1)


##### BUILD SCATTERPLOTS

wgbs.tissue = plot_scatter_pearsonCor_annotated(ge.prom, 'meth.diff', 'ge.log2FoldChange', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)

mbd.tissue = plot_scatter_pearsonCor_annotated(ge.prom, 'mbd.log2FoldChange', 'ge.log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)


mr.tissue = plot_scatter_pearsonCor_annotated(ge.prom, 'mr.log2FoldChange', 'ge.log2FoldChange', 'mdRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lty=2)



scatter_list = list(wgbs.tissue, mbd.tissue, mr.tissue)
mod_scatters = lapply(scatter_list, function(x) return(x + theme(axis.title = element_blank())))
scatter_plts = plot_grid(plotlist = mod_scatters, nrow=1)
vylab = ggdraw() + draw_label(bquote('-log'[10]~'p-value'), angle=90)
vplts = plot_grid(vylab, volc_plts, rel_widths = c(1,12))
sylab = ggdraw() + draw_label('GE difference', angle=90)
splts = plot_grid(sylab, scatter_plts, rel_widths = c(1,12))
xlab = ggdraw() + draw_label('promoter methylation difference')
plot_grid(vplts, splts, xlab, nrow = 3, rel_heights = c(12,12,1))


#plot only 


# SEESAWS -----------------------------------------------------------------

#replace the WGBS data with the exon data in gbm.dat as above
gbm.dat2 = eg.dat %>% 
  dplyr::select(name, l.fracMeth) %>% 
  right_join(gbm.dat, by = 'name') %>% 
  mutate(l.fracMeth = l.fracMeth.x)

#re-load response
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/tissue_promoter_response.Rdata')
ll  
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/genotype_promoter_response.Rdata')
ll

########## merge with genotype differences ##########
gbm_lvl_change = g.dat %>% 
  left_join(gbm.dat2, by = 'name') %>% 
  dplyr::select(name, l.fracMeth, meth.diff, mbd.score, mbd.log2FoldChange, mrB, mr.log2FoldChange) %>% 
  as_tibble()
gbm_lvl_ge_chnage = g.res %>% 
  left_join(gbm.dat2, by = 'name') %>% 
  dplyr::select(name, l.fracMeth, mbd.score, mrB, log2FoldChange) %>% 
  as_tibble()
########## merge with tissue differences ##########
gbm_lvl_change = t.dat %>% 
  left_join(gbm.dat2, by = 'name') %>% 
  dplyr::select(name, l.fracMeth, meth.diff, mbd.score, mbd.log2FoldChange, mrB, mr.log2FoldChange) %>% 
  as_tibble()
gbm_lvl_ge_chnage = t.res %>% 
  left_join(gbm.dat2, by = 'name') %>% 
  dplyr::select(name, l.fracMeth, mbd.score, mrB, log2FoldChange) %>% 
  as_tibble()
wgbs_xlims = c(log(0.01, 2), log(1, 2))
wgbs_ylims = c(-3, 3)
mbd_xlims = c(-2, 3)
mbd_ylims = c(-1,1)
mr_xlims = c(-4, 6)
mr_ylims = c(-3,2)
ge_ylims = c(-0.75, 0.75)
##################################################

seesaw_alpha = 0.02

#-------- GBM CHNAGE VS GBM LVL PLOTS
#WGBS
wgbs_ss = plot_scatter_pearsonCor_annotated(gbm_lvl_change, 'l.fracMeth', 'meth.diff', '', '', 
                                            ALPHA=seesaw_alpha, 
                                            xlim = wgbs_xlims,
                                            ylim = wgbs_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank()) +
  labs(title='WGBS') +
  scale_x_continuous(labels = log2_to_percent,
                     limits = wgbs_xlims)

#MBD-seq
mbd_ss = plot_scatter_pearsonCor_annotated(gbm_lvl_change, 'mbd.score', 'mbd.log2FoldChange', '', '',
                                           ALPHA=seesaw_alpha,
                                           xlim = mbd_xlims,
                                           ylim = mbd_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank()) +
  labs(title='MBD')

#mdRAD
mr_ss = plot_scatter_pearsonCor_annotated(gbm_lvl_change, 'mrB', 'mr.log2FoldChange', '', '',
                                          ALPHA=seesaw_alpha,
                                          xlim = mr_xlims,
                                          ylim = mr_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank()) +
  labs(title='mdRAD') 

#-------- GE CHNAGE VS GBM LVL PLOTS
#WGBS
wgbs_ss_ge = plot_scatter_pearsonCor_annotated(gbm_lvl_ge_chnage, 'l.fracMeth', 'log2FoldChange', '', '',
                                               ALPHA=seesaw_alpha,
                                               xlim = wgbs_xlims,
                                               ylim = ge_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  scale_x_continuous(labels = log2_to_percent) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank()) +
  scale_x_continuous(labels = log2_to_percent,
                     limits = wgbs_xlims)

#MBD-seq
mbd_ss_ge = plot_scatter_pearsonCor_annotated(gbm_lvl_ge_chnage, 'mbd.score', 'log2FoldChange', '', '',
                                              ALPHA=seesaw_alpha,
                                              xlim = mbd_xlims,
                                              ylim = ge_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank())

#mdRAD
mr_ss_ge = plot_scatter_pearsonCor_annotated(gbm_lvl_ge_chnage, 'mrB', 'log2FoldChange', '', '',
                                             ALPHA=seesaw_alpha,
                                             xlim = mr_xlims,
                                             ylim = ge_ylims) +
  geom_hline(yintercept = 0, lty=2) +
  geom_smooth(method='lm', se=FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  theme(axis.title = element_blank())

#------- ASSEMBLE
top_pans = plot_grid(wgbs_ss, mbd_ss, mr_ss, nrow=1)
bottom_pans = plot_grid(wgbs_ss_ge, mbd_ss_ge, mr_ss_ge, nrow=1)
topy = ggdraw() + draw_label('gbM difference', angle=90)
bottomy = ggdraw() + draw_label('GE difference', angle=90)
xlab =ggdraw() + draw_label('gbM level')
top = plot_grid(topy, top_pans, nrow=1, rel_widths = c(1,12))
bottom = plot_grid(bottomy, bottom_pans, nrow=1, rel_widths = c(1,12))
final = plot_grid(top, bottom, xlab, nrow=3, rel_heights = c(12,12,1))
final


# SEESAW DENSITY PLOTS ----------------------------------------------------


gbm_lvl_change

plot_amil_seesaw_density = function(dat, lvl_col, change_col, cut, xlab){
  dat = data.frame(dat)
  dat$meth_class = if_else(dat[,lvl_col] > cut,
                           'methylated',
                           'unmethylated')
  dat$meth_class = factor(dat$meth_class, levels=c('unmethylated','methylated'))
  plt=dat %>% 
    filter(!is.na(meth_class)) %>% 
    ggplot(aes_string(x=change_col)) +
    geom_density() +
    geom_vline(xintercept = 0, lty=2) +
    labs(y='', x=xlab) +
    facet_wrap(~meth_class, scales='free')
  return(plt)
}

#gbm
wgbs_dens = plot_amil_seesaw_density(dat=gbm_lvl_change, lvl_col='l.fracMeth', change_col='meth.diff', cut=-4, xlab='WGBS')
mbd_dens = plot_amil_seesaw_density(gbm_lvl_change, 'mbd.score', 'mbd.log2FoldChange', -0.25, 'MBD')
mr_dens = plot_amil_seesaw_density(gbm_lvl_change, 'mrB', 'mr.log2FoldChange', 0.4, 'mdRAD')
plot_grid(wgbs_dens, mbd_dens, mr_dens, nrow=1)

#ge
wgbs_dens = plot_amil_seesaw_density(dat=gbm_lvl_ge_chnage, lvl_col='l.fracMeth', change_col='log2FoldChange', cut=-4, xlab='WGBS')
mbd_dens = plot_amil_seesaw_density(gbm_lvl_ge_chnage, 'mbd.score', 'log2FoldChange', -0.25, 'MBD')
mr_dens = plot_amil_seesaw_density(gbm_lvl_ge_chnage, 'mrB', 'log2FoldChange', 0.4, 'mdRAD')
plot_grid(wgbs_dens, mbd_dens, mr_dens, nrow=1)




# things below here are kept only for reference ---------------------------


# FOR 500 BP WINDOWDS -----------------------------------------------------

ll=load('comparisons/datasets/tissue_500bp_response.Rdata')
ll  
ll=load('comparisons/datasets/genotype_500bp_response.Rdata')
ll  


#GET CLOSEST GENE


gdat = read_tsv('metadata/geneBounds.tsv',
                col_names = c('chr', 'start', 'end', 'description'))

CUT=0.2
bdat = t.dat %>% 
  dplyr::select(chr:end, meth.diff, pvalue, mr.log2FoldChange, mr.pvalue, mbd.log2FoldChange, mbd.pvalue) %>% 
  as_tibble() %>% 
  filter(pvalue < CUT,
         mr.pvalue < CUT,
         mbd.pvalue < CUT)
bdat

#get unique chrs
gchrs = unique(gdat$chr)
bchrs = unique(bdat$chr)
uchrs = gchrs[gchrs %in% bchrs]

for (chrom in uchrs){
  print(paste(chrom, '..', sep='.'))
  gsub = gdat %>% 
    filter(chr == chrom)
  bsub = bdat %>% 
    filter(chr == chrom)
  rdf = data.frame()
  
  for (i in 1:nrow(bsub)){
    s = as.numeric(bsub[i,'start'])
    e = as.numeric(bsub[i,'end'])
    leftGenes = gsub %>% 
      filter(start <= s) %>% 
      mutate(dist=start - s)
    rightGenes = gsub %>% 
      filter(start >= e) %>% 
      mutate(dist = start + e)
    closeLeft = leftGenes %>% 
      filter(dist == min(dist))
    closeRight = rightGenes %>% 
      filter(dist == min(dist))
    rres =bsub[i,]
    if (nrow(closeLeft)>0){
      rres$closeLeft = closeLeft$description
      rres$leftDist = closeLeft$dist 
    } else{
      rres$closeLeft = 'none'
      rres$leftDist = 'none'
    }
    if (nrow(closeRight)>0){
      rres$closeRight = closeRight$description
      rres$rightDist = closeRight$dist
    } else{
      rres$closeRight = 'none'
      rres$rightDist = 'none'
    }
    rdf = rbind(rdf, rres)
  }
}

split_out_name = function(sVector){
  s1 = sapply(sVector, function(x) strsplit(x, ';')[[1]][1])
  r = sub('ID=', '', s1)
  return(r)
}

rdf$name = split_out_name(rdf$closeRight)  
# rdf$name = split_out_name(rdf$closeLeft)

rdat = rdf %>% 
  left_join(t.res, by= 'name') %>% 
  dplyr::select(chr:mbd.pvalue, name:padj)


absdat = rdat %>% 
  mutate(ge = abs(log2FoldChange),
         pm = abs(meth.diff),
         mbd = abs(mbd.log2FoldChange),
         mr = abs(mr.log2FoldChange)) %>% 
  dplyr::select(ge:mr)

absdat %>% 
  ggplot(aes(x=pm, y=ge)) +
  geom_point(alpha=0.2)

wgbs.tissue = plot_scatter_pearsonCor_annotated(absdat, 'pm', 'ge', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

mbd.tissue = plot_scatter_pearsonCor_annotated(ge.prom, 'mbd.log2FoldChange', 'log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

mr.tissue = plot_scatter_pearsonCor_annotated(ge.prom, 'mr.log2FoldChange', 'log2FoldChange', 'MethylRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MethylRAD') +
  theme(plot.title = element_text(hjust=0.5))





#amil_closest_gene.R
source('invert_meth_ge_functions.R')

# LOAD EXPRESSION DATA --------------------------------------------------

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


# LOAD METHYLATION RESPONSE MEASURES FROM BENCHMARKING --------------------------------------
ll = load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/genotype_gbm_response.Rdata')
head(g.gbm.dat)
ll = load('bioprojects/amillepora_PRJNA601565/meth_response/other_methods/tissue_gbm_response.Rdata')
head(t.gbm.dat)




# LOAD NEW PICOMETHYL DIFFERENCE MEASURES -----------------------------------------------------
ll=load('bioprojects/amillepora_PRJNA601565/meth_response/genotype_genes_methylKit.Rdata')
ll
g.gbm = data.frame(myDiff) %>% 
  left_join(wbounds) %>% 
  left_join(g.res, by = 'name') %>% 
  as_tibble()

ll=load('bioprojects/amillepora_PRJNA601565/meth_response/tissue_genes_methylKit.Rdata')
ll
t.gbm = data.frame(myDiff) %>% 
  left_join(wbounds) %>% 
  left_join(t.res, by = 'name') %>% 
  as_tibble()




# COMPARE -----------------------------------------------------------------

#overall gbm
plot_scatter_pearsonCor_annotated(g.gbm, 'meth.diff', 'log2FoldChange', 'GBM change', 'GE change', ALPHA=0.1) +
  geom_smooth(method='lm', se=FALSE)






# LOAD CLOSEST GENE DIFFERNEITAL METHYLATION -------------------------------------------

read_closest_gene = function(filePath, maxDist, geneId ='ID='){
  mres = read_tsv(filePath)
  chrMatch = mres$gene.chr==mres$chr
  mres2 = mres %>% 
    filter(gene.distance <= maxDist) %>% 
    dplyr::select(-gene.chr, -pvalue)
  split1 = sapply(mres2$gene.description, function(x) strsplit(x, geneId, fixed=TRUE)[[1]][2])
  geneName = sapply(split1, function(x) strsplit(x, ';', fixed=TRUE)[[1]][1])
  mres2 %>% 
    mutate(name=geneName) %>% 
    dplyr::select(-chr, -start, -end, -strand, -gene.description, -gene.withinWindow, -gene.direction)
}

#function to merge methylkit and deseq results
merge_methylkit_deseq = function(mres, res){
  mdf = mres %>% 
    left_join(res, by='name')
}

#read in and merge
MAX_DIST=50e3
g.mres = 
  read_closest_gene('bioprojects/amillepora_PRJNA601565/meth_response/significant/genotype_500bpWindows_significant_closeGeneDistances.tsv',
                           MAX_DIST) %>% 
  merge_methylkit_deseq(res=g.res)
  
t.mres = read_closest_gene('bioprojects/amillepora_PRJNA601565/meth_response/significant/tissue_500bpWindows_significant_closeGeneDistances.tsv',
                  MAX_DIST) %>% 
  merge_methylkit_deseq(res=t.res)



#look at relationships
mres=g.mres
# mres=t.mres

MAX_DIST=2e3
mres %>% 
  filter(gene.distance < MAX_DIST) %>% 
  ggplot(aes(x=meth.diff, y=log2FoldChange)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method='lm')


#difference between those that are close and those that arnt

closeGenes = mres %>% 
  filter(gene.distance < MAX_DIST,
         !is.na(meth.diff)) %>% 
  pull(name)
length(closeGenes)

g.res %>% 
  mutate(isClose = name %in% closeGenes) %>% 
  ggplot(aes(x=isClose, y=abs(log2FoldChange))) +
  geom_boxplot() +
  lims(y=c(0,0.5))



#plot for series of distances from window to gene start

ll=load('~/gitreps/benchmarking_coral_methylation/comparisons/datasets/gbmLvl.Rdata')
head(gbm.dat)
keep = gbm.dat$name[gbm.dat$fracMeth < 0.04]


#filter for only closest genes
mres=g.mres %>% 
  filter(gene.enclosesWindow)


#filter for only methylated and closest genes
mres=g.mres %>% 
  filter(name %in% keep,
         gene.closest==TRUE,
         !gene.enclosesWindow)
  


dists = seq(1000, 10000, by=1000)
boxlist = list()
scatlist = list()
for (d in dists){
  closeGenes = mres %>% 
    filter(gene.distance < d,
           !is.na(meth.diff)) %>% 
    pull(name)
  
  box=mres %>% 
    mutate(isClose = name %in% closeGenes,
           abslog2FoldChange = abs(log2FoldChange)) %>% 
    ggplot(aes(x=isClose, y=abslog2FoldChange)) +
    geom_boxplot() +
    labs(subtitle = d) +
    lims(y=c(0, 0.45))
  
  scatter = mres %>% 
    filter(name %in% closeGenes) %>% 
    ggplot(aes(x=meth.diff, y=log2FoldChange)) +
    geom_point() +
    geom_smooth(method='lm') +
    labs(subtitle = d)
  
  
  print('----')
  print(d)
  print(length(closeGenes))
  boxlist[[as.character(d)]]=box
  scatlist[[as.character(d)]]=scatter
}

plot_grid(plotlist = boxlist)
plot_grid(plotlist = scatlist)










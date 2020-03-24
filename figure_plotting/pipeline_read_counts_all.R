#pipeline_read_counts_all.R

source('invert_meth_ge_functions.R')
bioprojectDirs


#function to read in pipeline counts
read_sdat = function(bpd){
  print(bpd)
  filePath = paste('bioprojects', bpd, 'rna_results/pipelineCounts.tsv', sep='/')
  sdat0 = read_tsv(filePath,
                   col_names=c('run', 'value', 'stat'))
  sdat = sdat0 %>% 
    filter(!grepl('_2', run)) %>% 
    mutate(run = sub('_dupsRemoved.bam', '', run),
           run = sub('.trim', '', run),
           run = sub('_1', '', run),
           run = sub('.', '-', run, fixed = TRUE),
           value = if_else(stat=='predupPropPaired',
                           value/2,
                           value),
           value = if_else(stat=='dedupPropPair',
                           value/2,
                           value),
           bioproject = bpd) %>% 
    mutate(stat = factor(stat, levels=c('rawCounts', 'trimmedCounts', 'predupPropPaired', 'dedupPropPair', 'geneCounted')))
  return(sdat)
}

sdatList=map(bioprojectDirs, read_sdat)
sdat = purrr::reduce(sdatList, rbind)

#plot barplot
bp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  ggplot(aes(x=stat, y=value, color=run, fill=run)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='Read count', x='Pipeline step', title='Moya experiment') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.5))


#plot scatter abs
lp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  ggplot(aes(x=stat, y=value, color=bioproject)) +
  geom_point() +
  geom_line(aes(group=run)) +
  theme(legend.position='none') +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot_grid(bp, lp)

  
#plot scatter prop raw
pp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  group_by(run) %>% 
  mutate(prop = value/max(value) ) %>% 
  ggplot(aes(x=stat, y=prop, color=bioproject)) +
    geom_point() +
    geom_line(aes(group=run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    # theme(legend.position='none') +
    lims(y=c(0,1))
l=cowplot::get_legend(pp)
plot_grid(lp, pp+theme(legend.position = 'none'), l, nrow=1, rel_widths = c(1,1,0.75))

cdat %>% 
  filter(!grepl('__', geneID))

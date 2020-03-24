
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
sdat0 = read_tsv('bioprojects/amillepora_PRJNA601565/rna_results/pipelineCounts/pipelineCounts.tsv',
                 col_names=c('run', 'value', 'stat'))
unique(sdat0$stat)
unique(sdat0$run)

sdat = sdat0 %>% 
  mutate(run = sub('.trim_dupsRemoved.bam', '', run),
         run = sub('.trim', '', run),
         run = sub('.', '-', run, fixed = TRUE)) %>% 
  mutate(stat = factor(stat, levels=c('rawCounts', 'trimmedCounts', 'predupMapped', 'dedupMapped', 'geneCounted')))

table(sdat$run)

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
  ggplot(aes(x=stat, y=value, color=run)) +
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
  ggplot(aes(x=stat, y=prop, color=run)) +
    geom_point() +
    geom_line(aes(group=run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none') +
    lims(y=c(0,1))

plot_grid(lp, pp)

cdat %>% 
  filter(!grepl('__', geneID))

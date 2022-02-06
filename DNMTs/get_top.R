library(tidyverse)



dat = read_tsv('./DNMTs/hits.tsv',
               col_names = c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                             'evalue', 'bitscore')) %>% 
  separate(saccver,
           into = c('species'),
           extra = 'drop',
           remove = FALSE)


top_hits = dat %>% 
  arrange(evalue) %>% 
  distinct(species, .keep_all = TRUE)
top_hits

top_hits %>% 
  select(saccver, pident, length, evalue) %>% 
  write_csv('DNMTs/top_hits.csv')

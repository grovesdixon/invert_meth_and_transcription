#!/usr/bin/env Rscript
#gene_length_from_gff.R


#This script is lazily done and should be improved for general use.


library(optparse)
option_list = list(
  
  make_option(c("--cds2mrna"), type="character", default=NULL, 
              help="Name of input gff"),

  make_option(c("--mrna2gene"), type="character", default='mRNA', 
              help="String indicating which feature you want parents for. Default = 'mRNA'"),

  make_option(c("--o"), type="character", default=NULL, 
              help="Name for output file")

)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cfile = opt$cds2mrna
mfile = opt$mrna2gene
out_name = opt$o



library(tidyverse)


c = read_tsv(cfile) %>%
  select(featureID, parentID)
m = read_tsv(mfile) %>%
  select(featureID, parentID)

colnames(c) = c('cds', 'rna')
colnames(m) = c('rna', 'gene')

c %>%
  left_join(m, by='rna') %>%
  select(-rna) %>%
  write_tsv(path=out_name)
#!/usr/bin/env Rscript
#gene_length_from_gff.R


#This script is lazily done and should be improved for general use.


library(optparse)
option_list = list(
  
  make_option(c("--gff"), type="character", default=NULL, 
              help="Name of input gff"),

  make_option(c("--feature"), type="character", default='mRNA', 
              help="String indicating which feature you want parents for. Default = 'mRNA'"),

  make_option(c("--idTag"), type="character", default='ID=', 
              help="String for splitting out the id tag. Default = 'ID='"),

  make_option(c("--parentTag"), type="character", default='Parent=', 
              help="String for splitting out the parent tag. Default = 'Parent='"),

  make_option(c("--o"), type="character", default=NULL, 
              help="Name for output file")

)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
gffIn = opt$gff
feature = opt$feature
idTag = opt$idTag
parentTag = opt$parentTag
outName = opt$o



library(tidyverse)
gdat = read.table(gffIn, sep="\t", quote="") %>% 
  filter(V3==feature) %>% 
  #get feature ID
  separate(V9, into=c('nothing1', 'fid0'), sep=idTag, remove=FALSE) %>% 
  separate(fid0, into=c('featureID'), sep=';') %>% 
  #get parentID
  separate(V9, into=c('nothing2', 'pid0'), sep=parentTag, remove=FALSE) %>% 
  separate(pid0, into=c('parentID'), sep=';') %>% 
  select(-nothing1, -nothing2) 

print("Gff head:")
print(head(gdat))


out = gdat[,c('V1', 'V3', 'V4', 'V5', 'V7', 'featureID', 'parentID')]
colnames(out) = c('chr', 'feature', 'start', 'stop', 'strand', 'featureID', 'parentID')
write.table(out, file=outName, sep="\t", quote=F, row.names=F)
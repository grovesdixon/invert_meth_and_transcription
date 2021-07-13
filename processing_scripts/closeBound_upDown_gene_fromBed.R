#!/usr/bin/env Rscript
#closest_upDown_gene_fromBed.R
library(tidyverse)
#script runs through lines in a bed file and returns genes with start sites
#(based on the first position in the 'gene' line in the gff with respect to the genes strand)
#It also returns the distance of the start site from the window start or window end
#depending on wether the gene was upstream or downstream of the window. If the starty site is
#within the window, distance is given with repsect to the window start


# PARSE ARGUMENTS ---------------------------------------------------------

suppressMessages(library(optparse))
option_list = list(
  
  make_option(c("--g"), type="character", default=NULL, 
              help="input gff"),
  
  make_option(c("--b"), type="character", default=NULL, 
              help="cutoff for raw pvalue (only those less will be output"),

  make_option(c("--c"), type="numeric", default=50000, 
              help="Cutoff for including a gene as 'close' to a window (default=50Kb)"),

    make_option(c("--header"), type="logical", default=TRUE, 
              help="logical for whethe rthe bedfile has a header. Defulat=TRUE"),
  
  make_option(c("--o"), type="character", default=NULL, 
              help="output prefix"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
bedIn = opt$b
gffIn = opt$g
outName = opt$o
CUT = opt$c
headerBool = opt$header


# bedIn = 'tissue_1KbWindows_significant.tsv'
# gffIn = 'Amil.coding.gff3'
# outName = 'testOut'
# CUT = 20000
# headerBool = TRUE



# READ IN DATA---------------------------------------------------------
print('Reading in gff...')
gdat = read_tsv(gffIn,
                col_names=c('chr','source','feature','start','end','x','strand','xx','description'),
                cols(chr=col_character(),
                     start=col_double(),
                     end=col_double(),
                     startSite=col_double())) %>% 
  filter(feature=='gene') %>% 
  mutate(startSite = if_else(strand=='+',
                             start,
                             end)) %>%
  select(chr,start,end,strand,startSite,description)

print('-------------')
print('gff header:')
print(head(gdat))

#read in the bedfile
print('Reading in bedfile...')
bdat = read.table(bedIn, stringsAsFactors=FALSE, header = headerBool)
#swap in first four column names, the rest can be whatever
colnames(bdat)[1:3] = c('chr', 'start', 'end')
bdat$start = as.numeric(bdat$start)
bdat$end = as.numeric(bdat$end)
print('-------------')
print('bed file header:')
print(head(bdat))


# ASSIGN CLOSE GENES -----------------------------------------------------
#get unique chrs
gchrs = unique(gdat$chr)
bchrs = unique(bdat$chr)
uchrs = gchrs[gchrs %in% bchrs]

#loop through unique chrs 
#for each window, get the closest upstream and downstream gene
print('Finding genes near windows...')
print(paste('cutoff =', CUT, 'bp'))
rdf = data.frame()
noGene.df = data.frame()
for (chrom in uchrs){
  print(paste(chrom, '..', sep='.'))
  gsub = gdat %>% 
    filter(chr == chrom)
  bsub = bdat %>% 
    filter(chr == chrom)
  #iterate through windows
  for (i in 1:nrow(bsub)){

    # print('------')
    # print(paste(i, ' of ', nrow(bsub), '...', sep=''))
    
    s = as.numeric(bsub[i,'start'])
    e = as.numeric(bsub[i,'end'])
    # bname = as.character(bsub[i,'name'])
    
    #subset genes located to left of window start or within the window
    #also assign distance from gene start to window start, and 
    #if gene is closest of all genes to left of the window
    gsub$withinWindow = gsub$startSite > s & gsub$startSite < e
    gsub$enclosesWindow = gsub$start < s & gsub$end > e
    toLeft = gsub %>% 
      filter(startSite <= s | withinWindow) %>% 
      mutate(distance=s - startSite,
             direction = 'leftOfWindowStart')
    minLeft = suppressWarnings(min(toLeft$distance, na.rm=TRUE))
    toLeft = toLeft %>% 
      mutate(closest = if_else(distance==minLeft,
             TRUE,
             FALSE))
    
    #repeat for genes located to right of the window end
    toRight = gsub %>% 
      filter(startSite >= e) %>% 
      mutate(distance = startSite - e,
             direction = 'rightOfWindowEnd')
    minRight = suppressWarnings(min(toRight$distance, na.rm=TRUE))
    toRight = toRight %>% 
      mutate(closest = if_else(distance==minRight,
             TRUE,
             FALSE))
    
    #bind the sets and filter for genes within the CUT threshold
    closeGenes = rbind(toLeft, toRight) %>% 
                 filter(distance < CUT)
    colnames(closeGenes) = paste('gene', colnames(closeGenes), sep='.')
    if (nrow(closeGenes) > 0){
      rres = closeGenes %>%
             mutate(chr=chrom,
                    start = s,
                    end = e) %>%
             select(append(c('chr','start','end'), colnames(closeGenes))) %>%
             data.frame()
      # print(head(rres))
      rdf = rbind(rdf, rres)
    } else{
      noGene.df = rbind(noGene.df, bsub[i,])
    }
  }
}

#now stick the results back onto the original bed file
print('------------')
print(head(rdf))
print('=============')
print(head(bdat))

finalDf = bdat %>%
              full_join(rdf, by=c('chr','start','end'))


#write out the results
finalDf %>% 
  write_tsv(path=outName)

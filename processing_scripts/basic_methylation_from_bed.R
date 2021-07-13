#!/usr/bin/env Rscript
#reduce_covSums.R

#PARSE ARUGMENTS
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(tidyverse))
suppressMessages(library(plotrix))
suppressMessages(library(EnvStats))

option_list = list(
  
  make_option(c("--cov"), type="character", default=NULL, 
              help="covInfile"),

  make_option(c("--bed"), type="character", default=NULL, 
              help="bedInfile"),

  make_option(c("--minCount"), type="integer", default=1, 
              help="Minimum number of reads a site must have to be counted"),

  make_option(c("--minRep"), type="double", default=0, 
              help="Minimum proportion of samples with filter passing data for the site to be kept"),

  make_option(c("--pFalsePos"), type="double", default=0.01, 
              help="Minimum proportion of samples with filter passing data for the site to be kept"),

  make_option(c("--methAlpha"), type="double", default=0.05, 
              help="Theshold for type 1 error for calling a site methylated given probability of false methylation call pFalsePos"),

  make_option(c("--o"), type="character", default='gbm_stats', 
              help="Output name")
)


print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
covFile = opt$cov
bedFile = opt$bed
minCount = opt$minCount
minRep = opt$minRep
propFalsePos = opt$pFalsePos
methAlpha = opt$methAlpha
outName = opt$o

#function to summarize stats
getStats = function(inputDF, gout){
        inputDF %>%
          summarize(chr=uchr,
            start = s,
            end = e,
            name = name,
            nM = sum(nM, na.rm=TRUE),
            nU = sum(nU, na.rm=TRUE),
            NcpgMeasured = n(),
            Nmcpg = sum(methylated),
            mCpG_per_CpG = Nmcpg/NcpgMeasured,
            fracMeth = sum(nM, na.rm=TRUE) / ( sum(nM, na.rm=TRUE) + sum(nU, na.rm=TRUE)),
            mnMeth = mean(pct.meth, na.rm=TRUE),
            medMeth = median(pct.meth, na.rm=TRUE),
            #geoMnMeth = geomMean(mPct, na.rm=TRUE),
            sdMeth = sd(pct.meth),
            stdErrMeth = std.error(pct.meth),
            maxMeth = max(pct.meth),
            minMeth = min(pct.meth))
      }


#READ IN DATA
print('Reading in bed file...')
bdat = read.table(bedFile, stringsAsFactors=FALSE)
colnames(bdat)=c('chr', 'start', 'end', 'name')
bdat = as_tibble(bdat)
print('Reading in cov file...')
cdat = read.table(covFile, stringsAsFactors=FALSE)
colnames(cdat)=c('chr', 'start', 'end', 'pct.meth', 'nM', 'nU', 'fileName')
cdat = as_tibble(cdat)
uchrs0 = unique(bdat$chr)
uchrs = uchrs0[uchrs0 %in% cdat$chr]
if (length(uchrs)==0){
  print('No gene regions found in this cov file.')
  print('Exiting')
  quit()
}

res = data.frame()
#LOOP THROUGH CHROMOSOMES AND GENES
for (chrNum in 1:length(uchrs)){
  uchr=uchrs[chrNum]
  print(paste(uchr, '...', sep=''))
  if (chrNum %% 100 == 0){
    print(paste('chr', chrNum, 'of', length(uchrs)))
  }
  csub = cdat %>% 
    filter(chr==uchr)
  bsub = bdat %>% 
    filter(chr==uchr)
  if (nrow(csub)==0){
    next
  }
  for (i in 1:nrow(bsub)){
    s=as.numeric(bsub[i,'start'])
    e=as.numeric(bsub[i,'end'])
    name=as.character(bsub[i,'name'])
    wsub = csub %>% 
      filter(start >= s,
             end <= e)
    if (nrow(wsub) > 0){
      dat=wsub

      #DO MIN COUNT FILTER
      # print("Filtering by read count...")
      f1 = dat %>% 
        mutate(tot=nM+nU) %>% 
        filter(tot>=minCount)

      before=nrow(dat)
      after=nrow(f1)
      pct = round(after/before, digits=3)*100
      #print(paste(c(pct, '% of sites passed minCount >= ', minCount), collapse=''))


      #DO MIN REP FILTER
      # print('Filtering by representation...')
      totSamples = length(unique(f1$fileName))
      keepSites = f1 %>% 
        group_by(chr, start, end) %>% 
        summarize(N=n(), rep=n()/totSamples, keep= (n()/totSamples) >=minRep) %>% 
        filter(keep)

      f2 = f1 %>% 
        filter(chr %in% keepSites$chr & start %in% keepSites$start)

      before = nrow(f1)
      after = nrow(f2)
      pct = round(after/before, digits=3)*100
      #print(paste(c(pct, '% of sites passed minRep >= ', minRep), collapse=''))


      #MAKE METHYLATION CALLS FOR EACH SAMPLE
      # print('Making site methylation calls...')

      f3 = f2 %>% 
        mutate(pFalsePos = unlist(map2(.x=nM,
                                       .y=tot,
                                       ~ binom.test(.x, .y, propFalsePos, alternative="greater")$p.value))
        )


      #CALL METH AND MERGE UP WITH GFF DATA FOR LENGTHS
      mdat = f3 %>% 
        mutate(methylated = pFalsePos < methAlpha,
               mPct = nM/(nM+nU))


      NcpgMeasured = length((mdat$start))
      Nmcpg = sum(mdat$methylated)

      #funciton to write out stats
      subres = getStats(mdat)
      res = rbind(res, subres)
    }
  }
}


print(paste('writing results to', outName))
write_tsv(res, path=outName)


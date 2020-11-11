#initialize_counts.R
#prepare reads for DESeq and get gist of data

library('readxl')
library('DESeq2')
source('invert_meth_ge_functions.R')
rm(list=ls())


rnaDir='bioprojects/primates/chimp/rna_results/'

#upload read counts from htseq
countsFile=paste(rnaDir, 'feature_counts_out.tsv', sep='')
counts = read.table(countsFile, header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]

#fix names
revisedNames = sub('_rna_dupsRemoved.bam', '', colnames(counts))
data.frame(colnames(counts), revisedNames)
colnames(counts) = revisedNames
head(counts)
length(unique(colnames(counts)))


#remove non-gene objects from counts data
dim(counts)
tots = apply(counts, 2, sum)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))


#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]



#SET UP COLDATA
sample = colnames(counts)


#read in data from sra
sra = read_excel('metadata/Blake_SraRunTable.xlsx') %>% 
  mutate(tissue = if_else(tissue=='heart (later reclassified as liver)',
                          'heartReclassedAsLiver',
                          tissue)) %>% 
  filter(`Assay Type`=='RNA-Seq',
         Organism == 'Pan troglodytes') %>% 
  dplyr::select(Individual, sex, tissue) %>% 
  mutate(Individual = paste("chimp", Individual, sep=''),
         Individual = sub('-', '.', Individual, fixed=TRUE),
         sample = paste(Individual, tissue, sep='_')) %>% 
  column_to_rownames('sample')

#set as coldata
coldata = sra[sample,]
coldata
dim(coldata)
sum(rownames(coldata)==colnames(counts))
write_tsv(coldata, path='bioprojects/primates/chimp/rna_results/rna_treat_table.txt')


#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))

#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds)

#get variance stabilized counts and save them
rld = rlog(dds)
rld.df=assay(rld)
colnames(rld.df) = colnames(counts)

#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have geneswith too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
rld=rld[rownames(rld.df),]
dim(rld.df)
dim(rld)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#build sample heatmaps 
library(pheatmap)
phm=pheatmap(cor(rld.df))

#plot pca
NTOP=round(nrow(counts)/10, digits=0)
pca_df = build_pca(rld.df, coldata, ntop = NTOP, pcs = 2)

addx=3
addy=2
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=tissue, shape = sex)) +
  geom_point(size=5) +
  lims(x=c(min(pca_df$PC1)-addx, max(pca_df$PC1+addx)),
       y=c(min(pca_df$PC2)-addy, max(pca_df$PC2+addy)))
  



#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#output
dim(counts)
dim(coldata)
deseqInfile=paste(rnaDir, 'deseqInput.Rdata', sep='')
rldFile=paste(rnaDir, 'rld.Rdata', sep='')
save(counts, coldata, file=deseqInfile)
save(rld.df, coldata, file=rldFile)


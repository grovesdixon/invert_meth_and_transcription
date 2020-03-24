#initialize_counts.R
#prepare reads for DESeq and get gist of data

rm(list=ls())
library('DESeq2')
source('invert_meth_ge_functions.R')



# SELECT THE BIOPROJECT ---------------------------------------------------

rnaDir='bioprojects/subcaste_bees_PRJNA104931/rna_results/'


# UPLOAD COUNTS -----------------------------------------------------------

#upload the counts
countsFile=paste(rnaDir, 'feature_counts_out.tsv', sep='')
counts = read.table(countsFile, header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]
revisedNames = sub('_dupsRemoved.bam', '', colnames(counts))
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

ttable=read.table('bioprojects/subcaste_bees_PRJNA104931/rna_treat_table.txt',
                  sep='\t',
                  header=TRUE, stringsAsFactors=FALSE)
rownames(ttable) = ttable$run
coldata=ttable[colnames(counts),]
remove(ttable)
coldata

#check
dim(coldata)
dim(counts)
sum(coldata$run==colnames(counts))==ncol(counts)

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
phm=pheatmap(cor(rld.df),
             labels_row = coldata$treatInfo)

#plot pca
NTOP=nrow(counts)
pca_df = build_pca(rld.df, coldata, ntop = NTOP, pcs = 2)
plot_rld_pca(pca_df,
             group_col = 'treatInfo',
             fix_coords = TRUE) + 
  theme(legend.position = 'right')
  

#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)



# SELECT 2 OUTLIER NURSES BASED ON PCA ------------------------------------

outlierNames = pca_df %>% 
  filter(treatInfo == 'Nurse',
         PC1 < 0) %>% 
  pull(run)
keepSamples = pca_df %>% 
  filter(!run %in% outlierNames) %>% 
  pull(run)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #number of samples left after outlier removal
print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
outlierNames
print(paste(nSamples, "samples were kept"))


#replot heatmap without outlier
rld.df = rld.df[, !colnames(rld.df) %in% outlierNames]
coldata = coldata[!rownames(coldata) %in% outlierNames, ]
pheatmap(cor(rld.df), labels_row =coldata$treatInfo)
plot_RNA_PCA(rld.df, coldata)




#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')
counts=counts[,!colnames(counts) %in% outlierNames]
coldata=coldata[!coldata$run %in% outlierNames,]
dim(counts)
dim(coldata)
deseqInfile=paste(rnaDir, 'deseqInput.Rdata', sep='')
rldFile=paste(rnaDir, 'rld.Rdata', sep='')
save(counts, coldata, file=deseqInfile)
save(rld.df, coldata, file=rldFile)


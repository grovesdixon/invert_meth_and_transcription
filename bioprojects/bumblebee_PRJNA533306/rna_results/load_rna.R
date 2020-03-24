#load_rna.R

rnaDir = 'bioprojects/bumblebee_PRJNA533306/rna_results/'

#upload the counts
countsFile=paste(rnaDir, 'feature_counts_out.tsv', sep='')
counts = read.table(countsFile, header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]
revisedNames = sub('_dupsRemoved.bam', '', colnames(counts))
colnames(counts) = revisedNames
head(counts)
length(unique(colnames(counts)))

#prepare treatment table
ttable=read.table('bioprojects/bumblebee_PRJNA533306/rna_treat_table.txt',
                  header=TRUE)
rownames(ttable) = ttable$run
coldata=ttable[colnames(counts),]
remove(ttable)
coldata


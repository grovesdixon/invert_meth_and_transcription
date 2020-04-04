#README

#overview
The repo includes 9 projects each from a different invertebrate species. Strings I use for these are:
"aiptasia", "stylophora", "silkworm", "termite", "maternal", "bumblebee", "subcaste", "roboust", "amillepora"

Data from each is in a subdirectory of bioprojects/
Each project includes a bisulfite dataset (meth) and an RNAseq dataset (GE)


#folder contents
assemble_counts_for_models.R -- R script used to assemble following list objects
ge_data_list.Rdata     -- list of dataframes with gene expression differences between treatment groups
meth_data_list.Rdata   -- list of dataframes with methylation counts for each gene in each sample with treatment labels
mtt_path_list.Rdata    -- list of treatment tables for the methylation analysis. (these are just for reference)
window_data_list.Rdata -- list of dataframes with methylation differences for significant 1Kb windows, with information of closeby genes


#details
ge_data_list.Rdata
	columns:
		name = gene name, should join with 'name' from other dataframes for that species
		log2FoldChange = treatment group 1 RNAseq counts over treatment group 0 on log2 scale (from DESEq2).
			(Treatment groups 1 and 0 match those given in meth_data_list)
		pvalue = p-value for the log2FoldChange (from DESeq2)

meth_data_list.Rdata
	columns:
		chr = the chromosome the gene is on
		start = gene's start position
		end = gene's end position
		name = gene name
		sample_num = arbitrary sample number
		numMeth = the number of methylated counts for the sample in the gene
		numUnmeth = the number of unmethylated counts for the sample in the gene
		treatment = binary indicator for sample's treatment group

mtt_list.Rdata
	columns:
		run = the cov file name (from Bismark)
		treatInfo = short description of treatment
		id = sample ID
		treat = binary indicator for sample's treatment group


window_data_list.Rdata
	This one is not very tidy. Because there were a LOT of 1kb windows,
	I only got the nearby gene information for the significant ones. So
	each of these dataframes only includes significant windows. 'Subcaste'
	didn't have significant windows so it's not included.
	columns:
		chr, start and end = (as in meth_data_list)
		pvalue = p-value for significance of the difference in methylation between treatment group 1 and 0 (from methylKit)(https://bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html)
		meth.diff = logit for methylation difference between group 1 and 0 (from methylKit)(https://bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html)
		gene.chr = the chromosome the gene is found on (should match chr for window)
		gene.start = gene's left boundary
		gene.end = gene's right boundary
		gene.strand = gene's strand
		gene.startSite = gene's start site (should be left boundary for + strand genes and right boundary for - strand genes)
		gene.withinWindow = boolean for whether the gene is enclosed in the window (unlikely for 1Kb windows)
		gene.enclosesWindow = boolean for whether the window falled entirely within the gene boundaries
		gene.distance = distance to the gene's start site
		gene.direction = the direction
		gene.closest = boolean for whether this is the closest gene to this window in that direction (leftOfWindowEnd or rightOfWindowEnd)
		name = gene's name (should join with name from other two dataframes)
		
	
	
	


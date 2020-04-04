#bioprojects/
This folder has the data files organized from each project.
The project directory names are given in bioprojects.txt.
Each subdirectory has:
bismark_summary -- data on how bismark ran
bisulfite_treat_table.txt -- table used to run methylKit
cdsID_geneID.tsv -- table relating CDS ids to gene ids from the .gff
meth_level -- tables of the methylation level for region types calcaulted by summing methylated and unmethylated counts from all samples across all CpGs within the region
meth_region_counts -- tables of filtered region counts generated with methylKit. These are the counts the differential methylation results are based on
meth_response -- results from methylKit getting differences in methylation between groups given in bisulfite_treat_table.txt
meth_response/significant/ -- only the significant 1Kb and 500bp windows (full files too big) along with closest genes (see closeBound_upDown_gene_fromBed.R)
rna_results/ -- all the RNAseq analysis files. Each project has its own copy of initialize_feature_counts.R which was run manually for QC
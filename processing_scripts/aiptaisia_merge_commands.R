

setwd("/Users/grovesdixon/lab_files/projects/mani_methylation/metadata/aiptasia_PRJNA415358")
x=read.table("aiptaisia_sample_mergin_input.yxy", header = T, stringsAsFactors=F)
head(x)

#DEAL WITH BISULFITE DATA
b=x[x$Assay_Type=="Bisulfite-Seq",]
nrow(b) #12 reps x 8 lanes

#assign biological replicate and lane
b$rep = sapply(b$Library_Name, function(x) strsplit(x, "_")[[1]][3])
b$lane = sapply(b$Library_Name, function(x) strsplit(x, "_")[[1]][4])
length(unique(b$rep))  #12 replicates
length(unique(b$lane)) #8 lanes each

#print catting commands to concatenate the bisulfite *.trim files
fileConn<-file('cat_bisulfite_trim_commands.txt')
reps = unique(b$rep)
all.commands = c()
for (r in reps){
	sub=b[b$rep==r,]
	runs = sub$Run
	fruns=paste(runs, "1.trim", sep="_")
	rruns=paste(runs, "2.trim", sep="_")
	forString = paste(fruns, collapse=" ")
	revString = paste(rruns, collapse=" ")
	forOut = paste(paste(">", r), "met_1.trim", sep='_')
	revOut = paste(paste(">", r), "met_2.trim", sep='_')
	commandFor = paste( paste("cat", forString), forOut)
	commandRev = paste( paste("cat", revString), revOut)
	print("---------")
	print(r)
	print(paste(r, "runs:"))
	print(runs)
	print("commands:")
	print(commandFor)
	print(commandRev)
	all.commands = append(all.commands, commandFor)
	all.commands = append(all.commands, commandRev)
}
writeLines(all.commands, fileConn)
close(fileConn)


#REPEAT FOR RNASEQ

rdf=x[x$Assay_Type=="RNA-Seq",]
nrow(rdf) #12 reps x 6 lanes

#assign biological replicate and lane
rdf$rep = sapply(rdf$Library_Name, function(x) strsplit(x, "_")[[1]][3])
rdf$lane = sapply(rdf$Library_Name, function(x) strsplit(x, "_")[[1]][5])
length(unique(rdf$rep))  #12 replicates
length(unique(rdf$lane)) #4 lanes used
table(rdf$lane)          #have some doubled up on lane 5 and 6


#print catting commands to concatenate the bisulfite *.trim files
fileConn<-file('cat_RNA_trim_commands.txt')
reps = unique(rdf$rep)
all.commands = c()
for (r in reps){
	sub=rdf[rdf$rep==r,]
	runs = sub$Run
	fruns=paste(runs, "1.trim", sep="_")
	rruns=paste(runs, "2.trim", sep="_")
	forString = paste(fruns, collapse=" ")
	revString = paste(rruns, collapse=" ")
	forOut = paste(paste(">", r), "rna_1.trim", sep='_')
	revOut = paste(paste(">", r), "rna_2.trim", sep='_')
	commandFor = paste( paste("cat", forString), forOut)
	commandRev = paste( paste("cat", revString), revOut)
	print("---------")
	print(r)
	print(paste(r, "runs:"))
	print(runs)
	print("commands:")
	print(commandFor)
	print(commandRev)
	all.commands = append(all.commands, commandFor)
	all.commands = append(all.commands, commandRev)
}
writeLines(all.commands, fileConn)
close(fileConn)









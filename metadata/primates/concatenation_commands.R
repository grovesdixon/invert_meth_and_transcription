#concatenation_commands.R
rm(list=ls())
library(tidyverse)
library(readxl)
setwd('~/projects/primate_gbm/')


sra = read_excel('SraRunTable.xlsx') %>% 
  mutate(tissue = if_else(tissue=='heart (later reclassified as liver)',
                          'heartReclassedAsLiver',
                          tissue)) ## just clean up the changed to liver thing


# get summary stats -------------------------------------------------------

#FOR RNA
rdat = sra %>% 
  filter(`Assay Type`=='RNA-Seq')
dim(rdat)

#4 runs per individual
rdat %>% 
  group_by(Individual, Organism) %>% 
  summarize(N=n())

#1 per tissue
rdat %>% 
  group_by(Individual, tissue) %>% 
  summarize(N=n()) %>% 
  pull(N) %>% 
  table()

#4 indvs * 3 speces * 4 tissues = 48

# OUTPUT TXT FILES FOR EACH ASSAY-ORGANISM  FOR RNA -------------------------------
unique(rdat$Organism)
unique(rdat$tissue)
select_cols = c('Run', 'Assay Type', 'Organism', 'sex', 'tissue')

#function to write subset
write_sub = function(df, organism, select_cols, out){
  df %>% 
    filter(Organism==organism) %>% 
    select(select_cols) %>% 
    write_tsv(out)
}

write_sub(rdat, "Homo sapiens", select_cols, 'metadata/human_rna.txt')
write_sub(rdat, "Pan troglodytes", select_cols, 'metadata/chimp_rna.txt')
write_sub(rdat, "Macaca mulatta", select_cols, 'metadata/macac_rna.txt')



# MAKE RENAMING COMMANDS FOR RNA ------------------------------------------

fileConn<-file('metadata/rename_rna_fastqs.txt')
reps = unique(rdat$Individual)
tissues = unique(rdat$tissue)
tissues
output_df = data.frame()
all.commands = c()
for (r in reps){
  sub0=rdat %>% 
    filter(Individual==r)
  for (t in tissues){
    sub = sub0 %>% 
      filter(tissue == t)
    runs = sub$Run
    org = unique(sub$Organism)
    sex = unique(sub$sex)
    ind_tissue_string = paste(r,t,sep='_')
    print(ind_tissue_string)
    print(runs)
    if (length(runs) < 1){
      next
    }
    fastqs = paste(runs, "1.fastq", sep="_")
    trims = paste(runs, "1.trim", sep="_")
    fqString = paste(fastqs, collapse=" ")
    trimString = paste(trims, collapse=" ")
    
    fqOut = paste(ind_tissue_string, "rna_1.fastq", sep='_')
    trimOut = paste(ind_tissue_string, "rna_1.trim", sep='_')
    commandFq = paste( paste("mv", fqString), fqOut)
    commandTrim = paste( paste("mv", trimString), trimOut)
    file_df = data.frame('ind_tissue' = ind_tissue_string,
                         'Assay Type' = 'WGBS',
                         'Organism' = org,
                         'sex' = sex,
                         'tissue' = t)
    print('file df:')
    print(file_df)
    output_df = rbind(output_df, file_df)
    print("---------")
    print(paste(r,t, sep='--'))
    print(paste(r, "runs:"))
    print(runs)
    print("commands:")
    print(commandFq)
    print(commandTrim)
    all.commands = append(all.commands, commandFq)
    all.commands = append(all.commands, commandTrim)
  }
}
writeLines(all.commands, fileConn)
close(fileConn)



# MAKE CONCATENATION COMMANDS FOR WGBS --------------------------------------

#FOR WGBS
bdat = sra %>% 
  filter(`Assay Type`=='Bisulfite-Seq')

bdat %>% 
  group_by(Individual) %>% 
  summarize(N=n())

#how many per tissue?
bdat %>% 
  group_by(Individual, tissue) %>% 
  summarize(N=n())


fileConn<-file('metadata/cat_bisulfite_fastqs.txt')
reps = unique(mod_bdat$Individual)
tissues = unique(mod_bdat$tissue)
tissues

output_df = data.frame()
all.commands = c()
sra = sra %>% 
  mutate(`Assay Type` = sub('-Seq', '', `Assay Type`))


for (r in reps){
  sub0=bdat %>% 
    filter(Individual==r)
  for (t in tissues){
    sub = sub0 %>% 
      filter(tissue == t)
    runs = sub$Run
    if (length(runs) < 1){
      next
    }
    org = unique(sub$Organism)
    sex = unique(sub$sex)
    ind_tissue_string = paste(r,t,sep='_')
    fastqs = paste(runs, "1.fastq", sep="_")
    trims = paste(runs, "1.trim", sep="_")
    fqString = paste(fastqs, collapse=" ")
    trimString = paste(trims, collapse=" ")
    fqOut = paste(paste(">", ind_tissue_string), "catted_1.fastq", sep='_')
    trimOut = paste(paste(">", ind_tissue_string), "catted_1.trim", sep='_')
    commandFq = paste( paste("cat", fqString), fqOut)
    commandTrim = paste( paste("cat", trimString), trimOut)
    file_df = data.frame('ind_tissue' = ind_tissue_string,
                         'Assay Type' = 'WGBS',
                         'Organism' = org,
                         'sex' = sex,
                         'tissue' = t)
    print('file df:')
    print(file_df)
    output_df = rbind(output_df, file_df)
  print("---------")
  print(paste(r,t, sep='--'))
  print(paste(r, "runs:"))
  print(runs)
  print("commands:")
  print(commandFq)
  print(commandTrim)
  all.commands = append(all.commands, commandFq)
  all.commands = append(all.commands, commandTrim)
  }
}
writeLines(all.commands, fileConn)
close(fileConn)


unique(sra$Organism)
#look at human individuals
sra %>% 
  filter(Organism=='Homo sapiens') %>% 
  pull(Individual) %>% 
  unique()

#look at chimp individuals
sra %>% 
  filter(Organism=="Pan troglodytes") %>% 
  pull(Individual) %>% 
  unique() %>% 
  data.frame() %>% 
  write_tsv('metadata/chimp_inds.txt', col_names = FALSE)

#look at chimp individuals
sra %>% 
  filter(Organism=="Macaca mulatta") %>% 
  pull(Individual) %>% 
  unique() %>% 
  data.frame() %>% 
  write_tsv('metadata/macac_inds.txt', col_names = FALSE)


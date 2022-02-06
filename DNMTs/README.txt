
#get the protein seqs
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/672/135/GCA_003672135.1_Obir_v5.4/GCA_003672135.1_Obir_v5.4_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/652/005/GCF_001652005.1_ASM165200v1/GCF_001652005.1_ASM165200v1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/696/155/GCF_000696155.1_ZooNev1.0/GCF_000696155.1_ZooNev1.0_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/905/235/GCF_014905235.1_Bmori_2016v1.0/GCF_014905235.1_Bmori_2016v1.0_protein.faa.gz
wget http://spis.reefgenomics.org/download/Spis.genome.annotation.pep.longest.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/417/965/GCA_001417965.1_Aiptasia_genome_1.1/GCA_001417965.1_Aiptasia_genome_1.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_protein.faa.gz


#make sure they are tagged
sed 's/>/>Aiptaisia_/' GCA_001417965.1_Aiptasia_genome_1.1_protein.faa > Aiptaisia.fa
sed 's/>/>Ant_/' GCA_003672135.1_Obir_v5.4_protein.faa > Ant.fa
sed 's/>/>Bumblebee_/' GCF_000214255.1_Bter_1.0_protein.faa > Bumblebee.fa
sed 's/>/>Termite_/' GCF_000696155.1_ZooNev1.0_protein.faa > Termite.fa
sed 's/>/>CarpenterBee_/' GCF_001652005.1_ASM165200v1_protein.faa > CarpenterBee.fa
sed 's/>/>HoneyBee_/' GCF_003254395.2_Amel_HAv3.1_protein.faa > HoneyBee.fa
sed 's/>/>Acropora_/' GCF_013753865.1_Amil_v2.1_protein.faa > Acropora.fa
sed 's/>/>Silkworm_/' GCF_014905235.1_Bmori_2016v1.0_protein.faa > Silkworm.fa
sed 's/>/>Stylophora_/' Spis.genome.annotation.pep.longest.fa > Stylophora.fa


#build concatenation and blast the human DNMT1 against it

#save blast results as hits.tsv

#select best hit for each species with get_top.R

#pull these out and align
#save as dnmt1s.aln

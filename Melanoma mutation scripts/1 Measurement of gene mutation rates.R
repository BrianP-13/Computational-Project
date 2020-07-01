library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: This script measures C>T mutation rates within protein coding genes


out<-"output_directory/"


## Retreive Ensembl gene id list with gneomic positions
EnsDb<-EnsDb.Hsapiens.v75
EnsDb_genes<-genes(EnsDb,columns = c(listColumns(EnsDb, "gene")))
EnsDb_genes$seq_coord_system<-NULL
EnsDb_genes$symbol<-NULL
seqlevelsStyle(EnsDb_genes)<-"UCSC"
EnsDb_genes<-keepSeqlevels(EnsDb_genes, value = standardChromosomes(EnsDb_genes),pruning.mode = "coarse")
genome(EnsDb_genes)<-"hg19"
EnsDb_genes<-sortSeqlevels(EnsDb_genes)
EnsDb_genes<-sort(EnsDb_genes)
EnsDb_genes
isSorted(EnsDb_genes)

# Keep genes with "protein_coding" bio_type descriptions
Prot_coding<-EnsDb_genes[which(EnsDb_genes$gene_biotype == "protein_coding"),] # Keep only protein coding genes
Prot_coding$gene_biotype<-NULL


## Somatic mutation counting

# Cancer mutation file
load(file = "file_path/MELA_AU_cutaneous_SSM.RData")
SSM<-SSM[,c(3:5,7,9)]

SSM<-makeGRangesFromDataFrame(SSM,keep.extra.columns = T)
seqlevelsStyle(SSM)<-"UCSC"
SSM<-sortSeqlevels(SSM)
SSM<-sort(SSM)
seqinfo(SSM)<-seqinfo(EnsDb_genes)


# C>T mutations
C_to_mutT<-SSM[which(SSM$reference_genome_allele=="C" ),]
C_to_mutT<-C_to_mutT[which(C_to_mutT$mutated_to_allele=="T"),]

# G>A mutations
G_to_mutA<-SSM[which(SSM$reference_genome_allele=="G"),]
G_to_mutA<-G_to_mutA[which(G_to_mutA$mutated_to_allele=="A"),]

## Calculate mutations per gene
Mut_Table_CtoT<-Prot_coding
Mut_Table_CtoT$C_to_T_count<-0

Mut_Table_GtoA<-Prot_coding
Mut_Table_GtoA$G_to_A_count<-0

for (i in 1:NROW(Prot_coding)){
  Gene_db<-Mut_Table_CtoT[i]
  hits<-findOverlaps(query = C_to_mutT,subject = Gene_db) # findOverlaps(query, subject)
  Overlaps<-as.data.frame(hits)
  Mut_Table_CtoT$C_to_T_count[i]<-nrow(Overlaps)
  
  Gene_db<-Mut_Table_GtoA[i]
  hits<-findOverlaps(query = G_to_mutA,subject = Gene_db) # findOverlaps(query, subject)
  Overlaps<-as.data.frame(hits)
  Mut_Table_GtoA$G_to_A_count[i]<-nrow(Overlaps)
}

## Calculate TOTAL C>T rate per gene
Mut_Table_CtoT$Normalized_C_to_T_rate<-0

Mut_Table_GtoA$Normalized_G_to_A_rate<-0

# Get Nucleotide content for hg19
BSgenome<-BSgenome.Hsapiens.UCSC.hg19
seqinfo(Mut_Table_CtoT)<-seqinfo(BSgenome)
seqlevels(Mut_Table_CtoT,pruning.mode="coarse")<-seqlevelsInUse(Mut_Table_CtoT)

DNA<-getSeq(x = BSgenome,names=Mut_Table_CtoT)
C_content<-letterFrequency(x = DNA, letters = "C")
G_content<-letterFrequency(x = DNA, letters = "G")

# Normalize mutation counts by nucleotide content and gene size (bp length)
Mut_Table_CtoT$Normalized_C_to_T_rate<-Mut_Table_CtoT$C_to_T_count/C_content
Mut_Table_GtoA$Normalized_G_to_A_rate<-Mut_Table_GtoA$G_to_A_count/G_content

# Combine all mut data 
mcols(Mut_Table_ALL)<-cbind(mcols(Mut_Table_ALL), mcols(Mut_Table_GtoA)[3:4])
Mut_Table_ALL$Total_C_to_T_rate<-Mut_Table_ALL$Normalized_C_to_T_rate+Mut_Table_ALL$Normalized_G_to_A_rate

# change column class to numeric type
Mut_Table_ALL$Normalized_C_to_T_rate<-as.numeric(Mut_Table_ALL$Normalized_C_to_T_rate)
Mut_Table_ALL$Normalized_G_to_A_rate<-as.numeric(Mut_Table_ALL$Normalized_G_to_A_rate)
Mut_Table_ALL$Total_C_to_T_rate<-as.numeric(Mut_Table_ALL$Total_C_to_T_rate)

save(Mut_Table_ALL, file = paste0(out, "Gene_mutation_table_ALL_CtoT_v2.RData"))



## End ##
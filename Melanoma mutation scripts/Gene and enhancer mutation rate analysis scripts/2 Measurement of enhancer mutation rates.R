library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: Measuring C>T mutation rates within IMR90 enhancer elements from the EnhancerAtlas database


out<-"output_directory/"


######
## Load IMR90 enhancer dataset
IMR90_Enh_pos<-import.bed(con = "file_path/IMR90_enhancer_list_EnhancerAtlas2.0.bed")
IMR90_Enh_pos<-as.data.frame(IMR90_Enh_pos)
IMR90_Enh_pos$start<-IMR90_Enh_pos$start-1
IMR90_Enh_pos<-makeGRangesFromDataFrame(IMR90_Enh_pos,keep.extra.columns = T)
IMR90_Enh_pos$name<-paste0("EnhAtlas-IMR90 enhancer ",seq(from=1,to=NROW(IMR90_Enh_pos)))
seqinfo(IMR90_Enh_pos)<-seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqlevels(IMR90_Enh_pos)<-seqlevelsInUse(IMR90_Enh_pos)


## Somatic mutation counting

# Cancer mutation file
load(file = "file_path/MELA_AU_cutaneous_SSM.RData")
SSM<-SSM[,c(3:5,7,9)]
SSM<-makeGRangesFromDataFrame(SSM,keep.extra.columns = T)
seqlevelsStyle(SSM)<-"UCSC"
SSM<-sortSeqlevels(SSM)
SSM<-sort(SSM)
SSM<-subset(SSM,seqnames!="chrM")
SSM<-subset(SSM,seqnames!="chrY")

seqlevels(SSM)<-seqlevelsInUse(SSM)
seqinfo(SSM)<-seqinfo(IMR90_Enh_pos)


## C>T mutations
C_to_mutT<-SSM[which(SSM$reference_genome_allele=="C" ),]
C_to_mutT<-C_to_mutT[which(C_to_mutT$mutated_to_allele=="T"),]

## G>A mutations
G_to_mutA<-SSM[which(SSM$reference_genome_allele=="G"),]
G_to_mutA<-G_to_mutA[which(G_to_mutA$mutated_to_allele=="A"),]


## Calculate mutations per gene
Mut_Table_CtoT<-IMR90_Enh_pos
Mut_Table_CtoT$C_to_T_count<-0

Mut_Table_GtoA<-IMR90_Enh_pos
Mut_Table_GtoA$G_to_A_count<-0

for (i in 1:NROW(IMR90_Enh_pos)){
  Enh_db<-Mut_Table_CtoT[i]
  hits<-findOverlaps(query = C_to_mutT,subject = Enh_db) # findOverlaps(query, subject)
  Overlaps<-as.data.frame(hits)
  Mut_Table_CtoT$C_to_T_count[i]<-nrow(Overlaps)
  
  Enh_db<-Mut_Table_GtoA[i]
  hits<-findOverlaps(query = G_to_mutA,subject = Enh_db) # findOverlaps(query, subject)
  Overlaps<-as.data.frame(hits)
  Mut_Table_GtoA$G_to_A_count[i]<-nrow(Overlaps)
}


## Calculate TOTAL C>T rate per gene
Mut_Table_CtoT$Normalized_C_to_T_rate<-0
Mut_Table_CtoT<-makeGRangesFromDataFrame(Mut_Table_CtoT,keep.extra.columns = T)

Mut_Table_GtoA$Normalized_G_to_A_rate<-0
Mut_Table_GtoA<-makeGRangesFromDataFrame(Mut_Table_GtoA,keep.extra.columns = T)


BSgenome<-BSgenome.Hsapiens.UCSC.hg19
seqinfo(Mut_Table_CtoT)<-seqinfo(BSgenome)
seqlevels(Mut_Table_CtoT,pruning.mode="coarse")<-seqlevelsInUse(Mut_Table_CtoT)

DNA<-getSeq(x = BSgenome,names=Mut_Table_CtoT)
C_content<-letterFrequency(x = DNA, letters = "C")
G_content<-letterFrequency(x = DNA, letters = "G")

# Normalize mutation counts by nucleotide content and gene size (bp length)
Mut_Table_CtoT$Normalized_C_to_T_rate<-Mut_Table_CtoT$C_to_T_count/C_content

Mut_Table_GtoA$Normalized_G_to_A_rate<-Mut_Table_GtoA$G_to_A_count/G_content

save(Mut_Table_CtoT, file = paste0(out, "IMR90_Enh_mutation_table_C_to_T.RData"))
save(Mut_Table_GtoA, file = paste0(out, "IMR90_Enh_mutation_table_G_to_A.RData"))


## Combine all mut data 
Mut_Table_ALL<-Mut_Table_CtoT

mcols(Mut_Table_ALL)<-cbind(mcols(Mut_Table_ALL), mcols(Mut_Table_GtoA)[2:3])
Mut_Table_ALL$Total_C_to_T_rate<-Mut_Table_ALL$Normalized_C_to_T_rate+Mut_Table_ALL$Normalized_G_to_A_rate

# change column class to numeric type
Mut_Table_ALL$Normalized_C_to_T_rate<-as.numeric(Mut_Table_ALL$Normalized_C_to_T_rate)
Mut_Table_ALL$Normalized_G_to_A_rate<-as.numeric(Mut_Table_ALL$Normalized_G_to_A_rate)
Mut_Table_ALL$Total_C_to_T_rate<-as.numeric(Mut_Table_ALL$Total_C_to_T_rate)

Enh_mutR<-as.data.frame(Mut_Table_ALL)

save(Enh_mutR, file = paste0(out, "IMR90_Enh_mutation_table_ALL_CtoT.RData"))



## END ##

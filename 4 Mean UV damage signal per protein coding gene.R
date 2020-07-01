library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)

# Description: This script will measure the average UV lesion signal (fold change) over every protien coding gene in the human genome.


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

isSorted(EnsDb_genes)

## Keep genes with "protein_coding" bio_type descriptions
Prot_coding<-EnsDb_genes[which(EnsDb_genes$gene_biotype == "protein_coding"),] # Keep only protein coding genes
Prot_coding$gene_biotype<-NULL


######
## Binning 6-4PP FC per gene

IP64bw<-"path_to_IP64_signal_track_file/IP64100vInput.fc.signal.bw"
IP64 <- import.bw(con= IP64_bw, as='RleList')
IP64[IP64 == 0] <- NA

seqinfo(IP64)<-seqinfo(Prot_coding)
seqlevels(Prot_coding)<-names(IP64)

IP64_gene<- binnedAverage(Prot_coding,numvar=IP64,varname='Mean_64PP_FC', na.rm = TRUE)

IP64_gene<-sortSeqlevels(IP64_gene)
IP64_gene<-sort(IP64_gene)

save(IP64_gene, file="output_directory/IP64_FC_per_prot_coding_gene.RData")



######
## Binning CPD FC per gene

CPD_bw<-'D:/Brian/Google Drive/Stanford University/LAB/Computational Project/Genomics Data/CPD100vInput.fc.signal.bw'
CPD <- import.bw(con= CPD_bw, as='RleList')
CPD[CPD == 0] <- NA

seqinfo(CPD)<-seqinfo(Prot_coding)
seqlevels(Prot_coding)<-names(CPD)

CPD_gene<- binnedAverage(Prot_coding,numvar=CPD,varname='Mean_CPD_FC', na.rm = TRUE)
CPD_gene<-sortSeqlevels(CPD_gene)
CPD_gene<-sort(CPD_gene)
isSorted(CPD_gene)

save(CPD_gene, file="output_directory/CPD_FC_per_prot_coding_gene.RData")


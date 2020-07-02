library(GenomicRanges)
#library(EnsDb.Hsapiens.v75)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(rtracklayer)

# Description: This script is for formatting the enhancer-target gene interaction dataset from EnhancerAtlas for IMR90. Will only keep protein coding genes


out<-"output_directory/"


## Enhancer-gene interactions
Enh_genes_RAW<-read.delim(file = "path_to_Enh_gene_dataset/EnhancerAtlas_IMR90_Enh_gene_interactions_v2.0.txt", header = F)

Enh_genes_df<-as.data.frame( matrix(nrow = nrow(Enh_genes_RAW),ncol = 9) )
colnames(Enh_genes_df)<-c("seqnames","start","end","Gene_ID","Gene_name","Gene_chr","Gene_TSS","Gene_strand","Confidence_score")

for (i in 1:nrow(Enh_genes_RAW)) {
  Split_col<-as.character(Enh_genes_RAW[i,1])
  
  Split_col<-unlist(strsplit(x = Split_col, split=":"))
  Split_col2<-unlist(strsplit(x = Split_col[2], split="-"))
  Split_col3<-unlist(strsplit(x = Split_col2[2], split="_"))
  Split_col4<-unlist(strsplit(x = Split_col3[2], split="[$]"))
  Split_col5<-unlist(strsplit(x = Split_col[2], split="[$]"))
  
  Enh_genes_df[i,1]<-Split_col[1] # seqnames
  Enh_genes_df[i,2]<-Split_col2[1] # start
  Enh_genes_df[i,3]<-Split_col3[1] # end
  Enh_genes_df[i,4]<-Split_col4[1] # Gene_ID
  Enh_genes_df[i,5]<-Split_col5[2] # Gene_name
  Enh_genes_df[i,6]<-Split_col5[3] # Gene_chr
  Enh_genes_df[i,7]<-Split_col5[4] # Gene_TSS
  Enh_genes_df[i,8]<-Split_col5[5] # Gene_strand
  Enh_genes_df[i,9]<-Enh_genes_RAW[i,2]
  print(paste0("Percent complete: ",signif(x = i/nrow(Enh_genes_RAW)*100,digits = 3),"%"))
}

## Keep only protein coding genes
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19)

EnsDb<-EnsDb.Hsapiens.v75
EnsDb_genes<-genes(EnsDb,columns = c(listColumns(EnsDb, "gene")))
EnsDb_genes$seq_coord_system<-NULL
EnsDb_genes$symbol<-NULL
seqlevelsStyle(EnsDb_genes)<-"UCSC"
EnsDb_genes<-keepSeqlevels(EnsDb_genes, value = standardChromosomes(EnsDb_genes),pruning.mode = "coarse")
genome(EnsDb_genes)<-"hg19"
EnsDb_genes<-sortSeqlevels(EnsDb_genes)

## Keep genes with "protein_coding" bio_type descriptions
Prot_coding<-EnsDb_genes[which(EnsDb_genes$gene_biotype == "protein_coding"),] # Keep only protein coding genes
Prot_coding$gene_biotype<-NULL

Gene_IDs<-unique(Prot_coding$gene_id)
IMR90_Enh_prot_genes<-Enh_genes_df

IMR90_Enh_prot_genes<-IMR90_Enh_prot_genes[which(IMR90_Enh_prot_genes$Gene_ID %in% Gene_IDs),]
rownames(IMR90_Enh_prot_genes)<-NULL

IMR90_Enh_prot_genes<-makeGRangesFromDataFrame(IMR90_Enh_prot_genes,keep.extra.columns = T)
save(IMR90_Enh_prot_genes,file = paste0(out, "IMR90_Enh_prot_coding_gene_interactions.RData"))



library(GenomicRanges)


## Description: Average together the total excision repair from BOTH strands ###


## Strategy:
# Step 1: Load data sets containing repair signal in protein coding genes (both plus and minus strand repair datasets)
# Step 2: Prepare forloop script to average repair rates from both strands.


out<-"output_directory/"


######
## 64PP repair time points

out_64<-out
XR64_path<-"path_to_merged_replicate_XRseq_files/"

file_names<-list.files(XR64_path) # 20 files
PatternList<-c("5min","20min","1h","2h","4h")


for (i in 1:5) {
  load(file = paste0(XR64_path, "XR64PP_PLUS_",  PatternList[i], "_Prot_genes.RData"))
  Pooled_plus_Prot<-as.data.frame(Pooled_plus_Prot)
  
  load(file = paste0(XR64_path, "XR64PP_MINUS_",  PatternList[i], "_Prot_genes.RData"))
  Pooled_minus_Prot<-as.data.frame(Pooled_minus_Prot)
  
  Mean_strand_repair<-Pooled_plus_Prot[,1:9]
  Mean_strand_repair$XR64_plus_mean<-Pooled_plus_Prot$XR64_plus_mean
  Mean_strand_repair$XR64_minus_mean<-Pooled_minus_Prot$XR64_minus_mean
  Mean_strand_repair$XR64_strand_mean<-rowMeans(Mean_strand_repair[,10:11])
  colnames(Mean_strand_repair)[12]<-paste0('XR64_strand_mean_', PatternList[i])
  
  Mean_strand_repair<-makeGRangesFromDataFrame(Mean_strand_repair, keep.extra.columns = T)
  save(Mean_strand_repair,file = paste0(out_64, "XR64PP_Strand_mean_",  PatternList[i], "_Prot_genes.RData"))
}


######
## CPD repair

out_CPD<-out
XRCPD_path<-"path_to_merged_replicate_XRseq_files/"

# File names
file_names<-list.files(XRCPD_path)
PatternList<-c("1h","4h","8h","16h","1d","2d")


for (i in 1:6) {
  load(file = paste0(XRCPD_path, "XRCPD_PLUS_",  PatternList[i], "_Prot_genes.RData"))
  Pooled_plus_Prot<-as.data.frame(Pooled_plus_Prot)
  
  load(file = paste0(XRCPD_path, "XRCPD_MINUS_",  PatternList[i], "_Prot_genes.RData"))
  Pooled_minus_Prot<-as.data.frame(Pooled_minus_Prot)
  
  Mean_strand_repair<-Pooled_plus_Prot[,1:9]
  Mean_strand_repair$XRCPD_plus_mean<-Pooled_plus_Prot$XRCPD_plus_mean
  Mean_strand_repair$XRCPD_minus_mean<-Pooled_minus_Prot$XRCPD_minus_mean
  Mean_strand_repair$XRCPD_strand_mean<-rowMeans(Mean_strand_repair[,10:11])
  colnames(Mean_strand_repair)[12]<-paste0('XRCPD_strand_mean_', PatternList[i])
  
  Mean_strand_repair<-makeGRangesFromDataFrame(Mean_strand_repair, keep.extra.columns = T)
  save(Mean_strand_repair,file = paste0(out_CPD, "XRCPD_Strand_mean_",  PatternList[i], "_Prot_genes.RData"))
}


## END ##
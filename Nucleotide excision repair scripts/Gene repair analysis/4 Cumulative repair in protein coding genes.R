library(GenomicRanges)


## Description: This script will generate the cummulative/total repair across all time points for each UV lesion. 
# There will be 2 types "Total" repair signals: Strnd average total repair, plus/minus strand total repair

## Strategy:
# First, I will calculate a total repair score using the averaged repair signal from both plus and minus strands
# Then I will calculate a total repair score for plus and minus strands respectively 


out<-"output_directory/"



### Total repair using BOTH plus and minus strands

## 64PP repair

out_64<-out
XR64_path<-"path_to_final_XRseq_files/"
file_names<-list.files(XR64_path) # 5 files
PatternList<-c("5min","20min","1h","2h","4h")

# Create basic database with protein coding genes BUT no repair data ( will add later )
load(file = paste0(XR64_path, "XR64PP_Strand_mean_",  PatternList[1], "_Prot_genes.RData"))
Mean_strand_repair<-as.data.frame(Mean_strand_repair)
Total_XR64_Prot_genes<-Mean_strand_repair[,1:9]

for (i in 1:5) {
  load(file = paste0(XR64_path, "XR64PP_Strand_mean_",  PatternList[i], "_Prot_genes.RData"))
  Mean_strand_repair<-as.data.frame(Mean_strand_repair)
  
  Total_XR64_Prot_genes<-cbind(Total_XR64_Prot_genes, Mean_strand_repair[,12])
  colnames(Total_XR64_Prot_genes)[9+i]<-paste0('XR64_strand_mean_', PatternList[i])
  
}

# Create a column for the sum total repair signal from all time points
Total_XR64_Prot_genes$Total_XR64_strand_mean_repair<-rowSums(x = Total_XR64_Prot_genes[,10:14], na.rm = T,dims = 1)
Total_XR64_Prot_genes<-makeGRangesFromDataFrame(Total_XR64_Prot_genes,keep.extra.columns = T)

save(Total_XR64_Prot_genes, file = paste0(out_64, 'Total_XR64_strand_mean_repair_Prot_genes.RData'))


## CPD repair

out_CPD<-out
XRCPD_path<-"path_to_final_XRseq_files/"
file_names<-list.files(XRCPD_path)
PatternList<-c("1h","4h","8h","16h","1d","2d")

# Create basic database with protein coding genes BUT no repair data ( will add later )
load(file = paste0(XRCPD_path, "XRCPD_Strand_mean_",  PatternList[1], "_Prot_genes.RData"))
Mean_strand_repair<-as.data.frame(Mean_strand_repair)
Total_XRCPD_Prot_genes<-Mean_strand_repair[,1:9]

for (i in 1:6) {
  load(file = paste0(XRCPD_path, "XRCPD_Strand_mean_",  PatternList[i], "_Prot_genes.RData"))
  Mean_strand_repair<-as.data.frame(Mean_strand_repair)
  
  Total_XRCPD_Prot_genes<-cbind(Total_XRCPD_Prot_genes, Mean_strand_repair[,12])
  colnames(Total_XRCPD_Prot_genes)[9+i]<-paste0('XRCPD_strand_mean_', PatternList[i])
  
}

# Create a column for the sum total repair signal from all time points
Total_XRCPD_Prot_genes$Total_XRCPD_strand_mean_repair<-rowSums(x = Total_XRCPD_Prot_genes[,10:14], na.rm = T,dims = 1)
Total_XRCPD_Prot_genes<-makeGRangesFromDataFrame(Total_XRCPD_Prot_genes,keep.extra.columns = T)

save(Total_XRCPD_Prot_genes, file = paste0(out_CPD, 'Total_XRCPD_strand_mean_repair_Prot_genes.RData'))


## END ##
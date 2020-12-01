library(GenomicRanges)


### Description: This script will generate the cummulative/total repair across all time points. 

## Strategy:
# First, I will calculate a total repair score using the averaged repair signal from both plus and minus strands

out<-"output_directory/"

######
### Total repair using BOTH plus and minus strands

## 64PP repair
path<-"path_to_final_XRseq_files/"
PatternList<-c("5min","20min","1h","2h","4h")

# Create basic database with chrom state regions
load(file = paste0(path, "XR64PP_Strand_mean_",  PatternList[1], "_core15.RData"))
Total_XR64_Chrom15<-as.data.frame(Mean_strand_repair)[,1:6]
head(Total_XR64_Chrom15)

for (i in 1:5) {
  load(file = paste0(path, "XR64PP_Strand_mean_",  PatternList[i], "_core15.RData"))
  Mean_strand_repair<-as.data.frame(Mean_strand_repair)
  
  Total_XR64_Chrom15<-cbind(Total_XR64_Chrom15, Mean_strand_repair[,9])
  colnames(Total_XR64_Chrom15)[6+i]<-paste0('XR64_strand_mean_', PatternList[i])
  
}
head(Total_XR64_Chrom15)

# Create a column for the sum total repair signal from all time points
Total_XR64_Chrom15$Total_XR64_strand_mean_repair<-rowSums(x = Total_XR64_Chrom15[,7:11], na.rm = T,dims = 1)
head(Total_XR64_Chrom15)

save(Total_XR64_Chrom15, file = paste0(out, 'Total_XR64_strand_mean_repair_core15.RData'))


## CPD repair
path<-"path_to_final_XRseq_files/"
PatternList<-c("1h","4h","8h","16h","1d","2d")

# Create basic database with chrom state regions
load(file = paste0(XRCPD_path, "XRCPD_Strand_mean_",  PatternList[1], "_core15.RData"))
Total_XRCPD_Chrom15<-as.data.frame(Mean_strand_repair)[,1:6]
head(Total_XRCPD_Chrom15)


for (i in 1:6) {
  load(file = paste0(XRCPD_path, "XRCPD_Strand_mean_",  PatternList[i], "_core15.RData"))
  Mean_strand_repair<-as.data.frame(Mean_strand_repair)
  
  Total_XRCPD_Chrom15<-cbind(Total_XRCPD_Chrom15, Mean_strand_repair[,9])
  colnames(Total_XRCPD_Chrom15)[6+i]<-paste0('XRCPD_strand_mean_', PatternList[i])
  
}
head(Total_XRCPD_Chrom15)

# Create a column for the sum total repair signal from all time points
Total_XRCPD_Chrom15$Total_XRCPD_strand_mean_repair<-rowSums(x = Total_XRCPD_Chrom15[,7:12], na.rm = T,dims = 1)

save(Total_XRCPD_Chrom15, file = paste0(out_CPD, 'Total_XRCPD_strand_mean_repair_core15.RData'))






### END ###

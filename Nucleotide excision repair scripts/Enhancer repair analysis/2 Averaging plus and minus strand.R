library(GenomicRanges)


## Description: Average together the total excision repair from BOTH strands

## Strategy:
# Step 1: Load data sets containing repair signal in protein coding genes (both plus and minus strand repair datasets)
# Step 2: Prepare forloop script to average repair rates from both strands.


out<-"output_directory/"


## Load IMR90 enhancer dataset
IMR90_Enh_pos<-import.bed(con = "file_path/IMR90_enhancer_list_EnhancerAtlas2.0.bed")
IMR90_Enh_pos<-as.data.frame(IMR90_Enh_pos)
IMR90_Enh_pos$start<-IMR90_Enh_pos$start-1
IMR90_Enh_pos<-makeGRangesFromDataFrame(IMR90_Enh_pos,keep.extra.columns = T)
IMR90_Enh_pos$name<-paste0("EnhAtlas-IMR90 enhancer ",seq(from=1,to=NROW(IMR90_Enh_pos)))

IMR90_Enh<-as.data.frame(IMR90_Enh_pos)


######
## 64PP repair time points

XR64_path<-"path_to_merged_replicate_XRseq_files/"
file_names<-list.files(XR64_path) # 20 files
PatternList<-c("5min","20min","1h","2h","4h")

Enh_XR64<-as.data.frame(IMR90_Enh_pos)

for (i in 1:5) {
  load(file = paste0(XR64_path, "XR64PP_PLUS_",  PatternList[i], "_IMR90_Enh.RData"))
  load(file = paste0(XR64_path, "XR64PP_MINUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
  Mean_strand_repair<-Pooled_plus_Enh[,1:6]
  Mean_strand_repair$XR64_plus_mean<-Pooled_plus_Enh$XR64_plus_mean
  Mean_strand_repair$XR64_minus_mean<-Pooled_minus_Enh$XR64_minus_mean
  Mean_strand_repair$XR64_strand_mean<-rowMeans(Mean_strand_repair[,7:8])
  Enh_XR64<-cbind(Enh_XR64,Mean_strand_repair[,9])
  colnames(Enh_XR64)[6+i]<-paste0('XR64_strand_mean_', PatternList[i])
  
}

# Measure total NER
Enh_XR64$Total_XR64_strand_mean_repair<-rowSums(x = Enh_XR64[,7:11], na.rm = T,dims = 1)

save(Enh_XR64, file = paste0(out,'Total_XR64_IMR90_Enh.RData'))


######
## CPD repair

XRCPD_path<-"path_to_merged_replicate_XRseq_files/"

# File names
file_names<-list.files(XRCPD_path)
PatternList<-c("1h","4h","8h","16h","1d","2d")

Enh_XRCPD<-as.data.frame(IMR90_Enh_pos)


for (i in 1:6) {
  load(file = paste0(XRCPD_path, "XRCPD_PLUS_",  PatternList[i], "_IMR90_Enh.RData"))
  load(file = paste0(XRCPD_path, "XRCPD_MINUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
  Mean_strand_repair<-Pooled_plus_Enh[,1:6]
  Mean_strand_repair$XRCPD_plus_mean<-Pooled_plus_Enh$XRCPD_plus_mean
  Mean_strand_repair$XRCPD_minus_mean<-Pooled_minus_Enh$XRCPD_minus_mean
  Mean_strand_repair$XRCPD_strand_mean<-rowMeans(Mean_strand_repair[,7:8])
  Enh_XRCPD<-cbind(Enh_XRCPD,Mean_strand_repair[,9])
  colnames(Enh_XRCPD)[6+i]<-paste0('XRCPD_strand_mean_', PatternList[i])
}

# Measure total NER
Enh_XRCPD$Total_XRCPD_strand_mean_repair<-rowSums(x = Enh_XRCPD[,7:12], na.rm = T,dims = 1)

save(Enh_XRCPD, file =paste0(out, 'Total_XRCPD_IMR90_Enh.RData'))


## END ##

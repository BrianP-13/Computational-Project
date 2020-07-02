library(GenomicRanges)
library(rtracklayer)


## Description: Measure NER rates within enhancer elements identified for IMR90 using only XRseq signals that uniquely map to the genome.


out<-"output_directory/"


## Load IMR90 enhancer dataset
IMR90_Enh_pos<-import.bed(con = "file_path/IMR90_enhancer_list_EnhancerAtlas2.0.bed")
IMR90_Enh_pos<-as.data.frame(IMR90_Enh_pos)
IMR90_Enh_pos$start<-IMR90_Enh_pos$start-1
IMR90_Enh_pos<-makeGRangesFromDataFrame(IMR90_Enh_pos,keep.extra.columns = T)
IMR90_Enh_pos$name<-paste0("EnhAtlas-IMR90 enhancer ",seq(from=1,to=NROW(IMR90_Enh_pos)))


######
## Step 2: Prepare for loop to process and bin XRseq data.

## XR64 repair data

out_64<-out
XR64_path<-"path_to_uniquely_mappable_XRseq_reads/"

# File names
file_names<-list.files(XR64_path) # 20 files
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
  IMR90_Enh<-IMR90_Enh_pos
  
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "Plus",value = T)
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "Plus",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "Minus",value = T)
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "Minus",value = T)
  
  ## Plus strand
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(XR64_path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep1_plus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep1_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_plus_Enh<-as.data.frame(Rep1_plus_Enh)
  head(Rep1_plus_Enh)
  print("Binning Rep1_plus_Enh complete")
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(XR64_path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep2_plus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  Rep2_plus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep2_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_plus_Enh<-as.data.frame(Rep2_plus_Enh)
  head(Rep2_plus_Enh)
  print("Binning Rep2_plus_Enh complete")
  
  # Combine both replicates and average them together
  Pooled_plus_Enh<-cbind(Rep1_plus_Enh, Rep2_plus_Enh[,7])
  Pooled_plus_Enh$R1_R2<-0
  Pooled_plus_Enh$R1_R2<-rowMeans(Pooled_plus_Enh[,c(7:8)], dims = 1)
  colnames(Pooled_plus_Enh)[c(7:9)]<-c("XR64_plus_Rep1","XR64_plus_Rep2","XR64_plus_mean")
  
  save(Pooled_plus_Enh,file = paste0(out_64, "XR64PP_PLUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(XR64_path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep1_minus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep1_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_minus_Enh<-as.data.frame(Rep1_minus_Enh)
  print("Binning Rep1_minus_Enh complete")
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(XR64_path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep2_minus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  Rep2_minus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep2_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_minus_Enh<-as.data.frame(Rep2_minus_Enh)
  print("Binning Rep2_minus_Enh complete")
  
  # Combine both replicates and average them together
  Pooled_minus_Enh<-cbind(Rep1_minus_Enh, Rep2_minus_Enh[,7])
  Pooled_minus_Enh$R1_R2<-0
  Pooled_minus_Enh$R1_R2<-rowMeans(Pooled_minus_Enh[,c(7:8)], dims = 1)
  colnames(Pooled_minus_Enh)[c(7:9)]<-c("XR64_minus_Rep1","XR64_minus_Rep2","XR64_minus_mean")
  
  save(Pooled_minus_Enh,file = paste0(out_64, "XR64PP_MINUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
}


######
## XR CPD repair data

out_CPD<-out
XRCPD_path<-"path_to_uniquely_mappable_XRseq_reads/"
list.files(XRCPD_path)

# File names
file_names<-list.files(XRCPD_path)
PatternList<-c("1h","4h","8h","16h","1d","2d")


for (i in 1:6) {
  IMR90_Enh<-IMR90_Enh_pos
  
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "Plus",value = T)
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "Plus",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "Minus",value = T)
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "Minus",value = T)
  
  ## Plus strand
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep1_plus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep1_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_plus_Enh<-as.data.frame(Rep1_plus_Enh)
  print("Binning Rep1_plus_Enh complete")
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep2_plus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  Rep2_plus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep2_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_plus_Enh<-as.data.frame(Rep2_plus_Enh)
  print("Binning Rep2_plus_Enh complete")
  
  # Combine both replicates and average them together
  Pooled_plus_Enh<-cbind(Rep1_plus_Enh, Rep2_plus_Enh[,7])
  Pooled_plus_Enh$R1_R2<-0
  Pooled_plus_Enh$R1_R2<-rowMeans(Pooled_plus_Enh[,c(7:8)], dims = 1)
  colnames(Pooled_plus_Enh)[c(7:9)]<-c("XRCPD_plus_Rep1","XRCPD_plus_Rep2","XRCPD_plus_mean")
  
  save(Pooled_plus_Enh,file = paste0(out_CPD, "XRCPD_PLUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep1_minus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep1_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_minus_Enh<-as.data.frame(Rep1_minus_Enh)
  print("Binning Rep1_minus_Enh complete")
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(IMR90_Enh)<-seqlevels(Rep2_minus_bw)
  IMR90_Enh<-sort(IMR90_Enh)
  
  Rep2_minus_Enh<-binnedAverage(bins = IMR90_Enh, numvar = Rep2_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_minus_Enh<-as.data.frame(Rep2_minus_Enh)
  print("Binning Rep2_minus_Enh complete")
  
  # Combine both replicates and average them together
  Pooled_minus_Enh<-cbind(Rep1_minus_Enh, Rep2_minus_Enh[,7])
  Pooled_minus_Enh$R1_R2<-0
  Pooled_minus_Enh$R1_R2<-rowMeans(Pooled_minus_Enh[,c(7:8)], dims = 1)
  colnames(Pooled_minus_Enh)[c(7:9)]<-c("XRCPD_minus_Rep1","XRCPD_minus_Rep2","XRCPD_minus_mean")
  
  save(Pooled_minus_Enh,file = paste0(out_CPD, "XRCPD_MINUS_",  PatternList[i], "_IMR90_Enh.RData"))
  
}


## END ##

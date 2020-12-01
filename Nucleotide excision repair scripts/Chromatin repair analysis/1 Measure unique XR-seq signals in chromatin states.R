library(GenomicRanges)
library(rtracklayer)

### Description: Measure NER rates within core 15 chromatin states for IMR90 using only XRseq signals that uniquely map to the genome ###

out<-"output_directory/"

# Load chromatin state file
core15<-file.path('file_path/E017_15_coreMarks_dense.bed.bgz')
GR_core15<-import.bed(core15)
mcols(GR_core15)$score<-NULL
mcols(GR_core15)$itemRgb<-NULL
mcols(GR_core15)$thick<-NULL


## XR64 repair data
out_64<-'D:/Brian/Google Drive/Stanford University/LAB/Computational Project/Genomics Data/XRseq processing/Uniquely mapped XRseq signal in chrom states for IMR90/NHF1_64/1 Merged replicates/'

# paths to 64PP XRseq bw files
path<-file.path("path_to_XRseq_files/")
#XR64_path<-'D:/Brian/Google Drive/Stanford University/LAB/Computational Project/Genomics Data/XRseq processing/Uniquely mapped XRseq signal in prot coding genes/NHF1_64/1 XR64 uniquely mapped regions/'

# File names
file_names<-list.files("path_to_XRseq_files/") # should be 20 files
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
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
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_plus_bw)
  GR_core15<-sort(GR_core15)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep1_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_plus_Chrom15<-as.data.frame(Rep1_plus_Chrom15)
  print(paste0("Binning Rep1_plus_Chrom15_", PatternList[i]," complete"))
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep2_plus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep2_plus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep2_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_plus_Chrom15<-as.data.frame(Rep2_plus_Chrom15)
  print(paste0("Binning Rep2_plus_Chrom15_", PatternList[i]," complete"))
  
  # Combine both replicates and average them together
  Pooled_plus_Chrom15<-cbind(Rep1_plus_Chrom15, Rep2_plus_Chrom15[,7])
  Pooled_plus_Chrom15$R1_R2<-0
  Pooled_plus_Chrom15$R1_R2<-rowMeans(Pooled_plus_Chrom15[,c(7:8)], dims = 1)
  colnames(Pooled_plus_Chrom15)[c(7:9)]<-c("XR64_plus_Rep1","XR64_plus_Rep2","XR64_plus_mean")
  
  Pooled_plus_Chrom15<-makeGRangesFromDataFrame(Pooled_plus_Chrom15, keep.extra.columns = T)
  Pooled_plus_Chrom15<-sortSeqlevels(Pooled_plus_Chrom15)
  Pooled_plus_Chrom15<-sort(Pooled_plus_Chrom15)
  
  save(Pooled_plus_Chrom15,file = paste0(out, "XR64PP_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_minus_bw)
  GR_core15<-sort(GR_core15)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep1_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_minus_Chrom15<-as.data.frame(Rep1_minus_Chrom15)
  print(paste0("Binning Rep1_minus_Chrom15_", PatternList[i]," complete"))
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep2_minus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep2_minus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep2_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_minus_Chrom15<-as.data.frame(Rep2_minus_Chrom15)
  print(paste0("Binning Rep2_minus_Chrom15_", PatternList[i]," complete"))
  
  # Combine both replicates and average them together
  Pooled_minus_Chrom15<-cbind(Rep1_minus_Chrom15, Rep2_minus_Chrom15[,7])
  Pooled_minus_Chrom15$R1_R2<-0
  Pooled_minus_Chrom15$R1_R2<-rowMeans(Pooled_minus_Chrom15[,c(7:8)], dims = 1)
  colnames(Pooled_minus_Chrom15)[c(7:9)]<-c("XR64_minus_Rep1","XR64_minus_Rep2","XR64_minus_mean")
  
  Pooled_minus_Chrom15<-makeGRangesFromDataFrame(Pooled_minus_Chrom15, keep.extra.columns = T)
  Pooled_minus_Chrom15<-sortSeqlevels(Pooled_minus_Chrom15)
  Pooled_minus_Chrom15<-sort(Pooled_minus_Chrom15)
  
  save(Pooled_minus_Chrom15,file = paste0(out, "XR64PP_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  print(paste0("ALL ", PatternList[i]," time points complete"))
}



## XR CPD repair data

# paths to CPD XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files(path)
PatternList<-c("1h","4h","8h","16h","1d","2d")


for (i in 1:6) {
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
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_plus_bw)
  GR_core15<-sort(GR_core15)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep1_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_plus_Chrom15<-as.data.frame(Rep1_plus_Chrom15)
  print(paste0("Binning Rep1_plus_Chrom15_", PatternList[i]," complete"))
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep2_plus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep2_plus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep2_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_plus_Chrom15<-as.data.frame(Rep2_plus_Chrom15)
  print(paste0("Binning Rep2_plus_Chrom15_", PatternList[i]," complete"))
  
  # Combine both replicates and average them together
  Pooled_plus_Chrom15<-cbind(Rep1_plus_Chrom15, Rep2_plus_Chrom15[,7])
  Pooled_plus_Chrom15$R1_R2<-0
  Pooled_plus_Chrom15$R1_R2<-rowMeans(Pooled_plus_Chrom15[,c(7:8)], dims = 1)
  colnames(Pooled_plus_Chrom15)[c(7:9)]<-c("XRCPD_plus_Rep1","XRCPD_plus_Rep2","XRCPD_plus_mean")
  
  Pooled_plus_Chrom15<-makeGRangesFromDataFrame(Pooled_plus_Chrom15, keep.extra.columns = T)
  Pooled_plus_Chrom15<-sortSeqlevels(Pooled_plus_Chrom15)
  Pooled_plus_Chrom15<-sort(Pooled_plus_Chrom15)
  
  save(Pooled_plus_Chrom15,file = paste0(out, "XRCPD_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_minus_bw)
  GR_core15<-sort(GR_core15)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep1_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_minus_Chrom15<-as.data.frame(Rep1_minus_Chrom15)
  print(paste0("Binning Rep1_minus_Chrom15_", PatternList[i]," complete"))
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep2_minus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep2_minus_Chrom15<-binnedAverage(bins = GR_core15, numvar = Rep2_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_minus_Chrom15<-as.data.frame(Rep2_minus_Chrom15)
  print(paste0("Binning Rep2_minus_Chrom15_", PatternList[i]," complete"))
  
  # Combine both replicates and average them together
  Pooled_minus_Chrom15<-cbind(Rep1_minus_Chrom15, Rep2_minus_Chrom15[,7])
  Pooled_minus_Chrom15$R1_R2<-0
  Pooled_minus_Chrom15$R1_R2<-rowMeans(Pooled_minus_Chrom15[,c(7:8)], dims = 1)
  colnames(Pooled_minus_Chrom15)[c(7:9)]<-c("XRCPD_minus_Rep1","XRCPD_minus_Rep2","XRCPD_minus_mean")
  
  Pooled_minus_Chrom15<-makeGRangesFromDataFrame(Pooled_minus_Chrom15, keep.extra.columns = T)
  Pooled_minus_Chrom15<-sortSeqlevels(Pooled_minus_Chrom15)
  Pooled_minus_Chrom15<-sort(Pooled_minus_Chrom15)
  
  save(Pooled_minus_Chrom15,file = paste0(out, "XRCPD_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  print(paste0("ALL ", PatternList[i]," time points complete"))
}


### END ###

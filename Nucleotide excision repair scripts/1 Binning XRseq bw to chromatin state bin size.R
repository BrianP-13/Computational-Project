library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: This script is the first step in measuring the excision repair levels for 6-4PP and CPD lesions within IMR90 chromatin states.
# This script also takes the additional step of averaging together both replicates (They are separate steps when measuring repair in 1Mb binned genome).


out<-"output_directory/"

### Binning XRseq repair signal to regions characterized by a single chromatin state ###
## Script Description: This R script will calculate the average repair signal for a particular 
## UV lesion type within regions of the genome characterized as one of the core 15 chromatin states

### Strategy: 
## Step 1: Bin ALL XR seq files to chromatin state bin sizes. Then combine the replicates for each
## time point and each strand by averaging the repair signal. Then save the averaged XRseq repair values
## Step 2: Combine the repair signals for both strands by averaging the signal from both strands (plus and minus)
## Step 3: Using the consolidated repair values previously calculated, create a composite plot overlapping
## all repair values for a specific time point and create a scatter plot with regression line.

##
# Load chromatin state file
core15<-file.path("file_path/E017_15_coreMarks_dense.bed.bgz")
GR_core15<-import.bed(core15)
mcols(GR_core15)$score<-NULL
mcols(GR_core15)$itemRgb<-NULL
mcols(GR_core15)$thick<-NULL


######
## Binning XRseq-64PP files to core 15 chromatinstate model

# paths to 64PP XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files("path_to_XRseq_files/") # should be 20 files
PatternList<-c("5min","20min","1h","2h","4h")

for (i in 1:5) {
  print(paste0("Working on ", PatternList[i], " XRseq files"))
  
  # Plus strand
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "PLUS",value = T)
  
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "PLUS",value = T)
  
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_plus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep1_plus_15<-binnedAverage(bins = GR_core15, numvar = Rep1_plus_bw, varname = "IP64 average", na.rm = TRUE)
  Rep1_plus_15<-as.data.frame(Rep1_plus_15)
  print("Binning Rep1_plus complete")
  
  Rep2_plus_15<-binnedAverage(bins = GR_core15, numvar = Rep2_plus_bw, varname = "IP64 average", na.rm = TRUE)
  Rep2_plus_15<-as.data.frame(Rep2_plus_15)
  print("Binning Rep2_plus complete")
  
  Pooled_plus<-cbind(Rep1_plus_15,Rep2_plus_15[,7])
  Pooled_plus$R1_R2<-0
  Pooled_plus$R1_R2<-rowMeans(Pooled_plus[,c(7:8)], dims = 1)
  colnames(Pooled_plus)[c(7:9)]<-c("Rep1","Rep2","mean_Rep")
  
  save(Pooled_plus,file = paste0(out, "XR64PP_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  
  rm(Rep1_plus)
  rm(Rep2_plus)
  rm(Rep1_plus_bw)
  rm(Rep2_plus_bw)
  rm(Rep1_plus_15)
  rm(Rep2_plus_15)
  rm(Pooled_plus)
  
  # Minus strand
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "MINUS",value = T)
  
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "MINUS",value = T)
  
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_minus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep1_minus_15<-binnedAverage(bins = GR_core15, numvar = Rep1_minus_bw, varname = "IP64 average", na.rm = TRUE)
  Rep1_minus_15<-as.data.frame(Rep1_minus_15)
  print("Binning Rep1_minus complete")
  
  Rep2_minus_15<-binnedAverage(bins = GR_core15, numvar = Rep2_minus_bw, varname = "IP64 average", na.rm = TRUE)
  Rep2_minus_15<-as.data.frame(Rep2_minus_15)
  print("Binning Rep2_minus complete")
  
  Pooled_minus<-cbind(Rep1_minus_15,Rep2_minus_15[,7])
  Pooled_minus$R1_R2<-0
  Pooled_minus$R1_R2<-rowMeans(Pooled_minus[,c(7:8)], dims = 1)
  colnames(Pooled_minus)[c(7:9)]<-c("Rep1","Rep2","mean_Rep")
  
  save(Pooled_minus,file = paste0(out, "XR64PP_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  
}


######
## Binning XRseq-CPD files to core 15 chromatinstate model

# paths to CPD XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files("path_to_XRseq_files/")
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:6) {
  print(paste0("Working on ", PatternList[i], " XRseq files"))
  
  # Plus strand
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "PLUS",value = T)
  
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "PLUS",value = T)
  
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_plus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep1_plus_15<-binnedAverage(bins = GR_core15, numvar = Rep1_plus_bw, varname = "CPD average", na.rm = TRUE)
  Rep1_plus_15<-as.data.frame(Rep1_plus_15)
  print("Binning Rep1_plus complete")
  
  Rep2_plus_15<-binnedAverage(bins = GR_core15, numvar = Rep2_plus_bw, varname = "CPD average", na.rm = TRUE)
  Rep2_plus_15<-as.data.frame(Rep2_plus_15)
  print("Binning Rep2_plus complete")
  
  Pooled_plus<-cbind(Rep1_plus_15,Rep2_plus_15[,7])
  Pooled_plus$R1_R2<-0
  Pooled_plus$R1_R2<-rowMeans(Pooled_plus[,c(7:8)], dims = 1)
  colnames(Pooled_plus)[c(7:9)]<-c("Rep1","Rep2","mean_Rep")
  
  save(Pooled_plus,file = paste0(out, "XRcpd_",  PatternList[i], "_R1_R2_PLUS_core15.RData"))
  
  rm(Rep1_plus)
  rm(Rep2_plus)
  rm(Rep1_plus_bw)
  rm(Rep2_plus_bw)
  rm(Rep1_plus_15)
  rm(Rep2_plus_15)
  rm(Pooled_plus)
  
  # Minus strand
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "MINUS",value = T)
  
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "MINUS",value = T)
  
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(GR_core15)<-seqlevels(Rep1_minus_bw)
  GR_core15<-sort(GR_core15)
  
  Rep1_minus_15<-binnedAverage(bins = GR_core15, numvar = Rep1_minus_bw, varname = "CPD average", na.rm = TRUE)
  Rep1_minus_15<-as.data.frame(Rep1_minus_15)
  print("Binning Rep1_minus complete")
  
  Rep2_minus_15<-binnedAverage(bins = GR_core15, numvar = Rep2_minus_bw, varname = "CPD average", na.rm = TRUE)
  Rep2_minus_15<-as.data.frame(Rep2_minus_15)
  print("Binning Rep2_minus complete")
  
  Pooled_minus<-cbind(Rep1_minus_15,Rep2_minus_15[,7])
  Pooled_minus$R1_R2<-0
  Pooled_minus$R1_R2<-rowMeans(Pooled_minus[,c(7:8)], dims = 1)
  colnames(Pooled_minus)[c(7:9)]<-c("Rep1","Rep2","mean_Rep")
  
  save(Pooled_minus,file = paste0(out, "XRcpd_",  PatternList[i], "_R1_R2_MINUS_core15.RData"))
  
}


## END ##
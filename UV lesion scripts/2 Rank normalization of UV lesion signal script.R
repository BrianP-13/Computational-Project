library(GenomicRanges)

## Description: This script will normalize the binned UV lesion signals to a normal distribution using ranked values and mapping them to a 
## normal distribution (Rank based Inverse Normal Transformation)

## CSV files containing 100kb and 1Mb binned genome with CPD and 6-4PP signal (columns 6 and 7):
# UVlesion_mastertable_100kb.csv
# UVlesion_mastertable_1Mb.csv

out<-"output_directory/"

# Rank normalization functions
rankitNormalize_vector <- function(x) {
  
  stopifnot(is.numeric(x))
  x <- qnorm((rank(x) - 0.5) / length(x))
  return(x)
}

rankitNormalize <- function(x, IND = 1) {
  
  # Normalizes rows (IND = 1) or columns (IND = 2)
  # to a quantile standard normalization (i.e. rankit)
  
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(x))
  
  rowNames <- rownames(x)
  colNames <- colnames(x)
  
  x <- apply(x, IND, rankitNormalize_vector)
  if(IND == 1)
    x <- t(x)
  
  rownames(x) <- rowNames
  colnames(x) <- colNames
  
  return(x)
  
}


######
## Rank normalization of 100kb binned signal

# Load 100kb binned UV lesion file and conduct normalization
Lesion_path<-"path_to_file_with_binned_64PP_and_CPD_signal/UVlesion_mastertable_100kb.csv"
UVLesion_MT<-read.csv(Lesion_path,header = T )
UVLesion_MT<-UVLesion_MT[,c(2:8)]
UVLesion_MT<-UVLesion_MT[which(UVLesion_MT$seqnames!="chrY"),]
UVLesion_MT<-UVLesion_MT[which(UVLesion_MT$seqnames!="chrM"),]
UVLesion_MT<-na.omit(UVLesion_MT)

CPD_64PP<-UVLesion_MT[,c(6,7)]
CPD_64PP<-as.matrix(CPD_64PP)

ranked_UV_lesions<-rankitNormalize(x= CPD_64PP, IND = 2)
ranked_UV_lesions<-as.data.frame(ranked_UV_lesions)

UVLesion_MT<-cbind(UVLesion_MT, ranked_UV_lesions)
colnames(UVLesion_MT)[c(8,9)]<-c("ranked_CPD","ranked_64PP")

# Convert dataframe to a GRanges type and sort according to chromosome order (e.g. chr1, chr2, ect)
UVLesion_MT<-makeGRangesFromDataFrame(UVLesion_MT, keep.extra.columns = T)
UVLesion_MT<-sortSeqlevels(UVLesion_MT)
UVLesion_MT<-sort(UVLesion_MT)
isSorted(UVLesion_MT)

# Measure the difference between ranked 6-4 and CPD
Ranked_diff<-UVLesion_MT
Ranked_diff$Difference<-Ranked_diff$ranked_64PP-Ranked_diff$ranked_CPD

# Change back to data.frame and save
Ranked_diff<-as.data.frame(Ranked_diff)
save(Ranked_diff, file = paste0(out, "UV lesions rank normalized to norm dist 100kb.RData"))


######
## Rank normalization of 1Mb binned signal

# Load 1Mb binned UV lesion file and conduct normalization
Lesion_path<-"path_to_file_with_binned_64PP_and_CPD_signal/UVlesion_mastertable_1Mb.csv"
UVLesion_MT<-read.csv(Lesion_path,header = T )
UVLesion_MT<-UVLesion_MT[,c(2:8)]
UVLesion_MT<-UVLesion_MT[which(UVLesion_MT$seqnames!="chrY"),]
UVLesion_MT<-UVLesion_MT[which(UVLesion_MT$seqnames!="chrM"),]
UVLesion_MT<-na.omit(UVLesion_MT)

CPD_64PP<-UVLesion_MT[,c(6,7)]
CPD_64PP<-as.matrix(CPD_64PP)

ranked_UV_lesions<-rankitNormalize(x= CPD_64PP, IND = 2)
ranked_UV_lesions<-as.data.frame(ranked_UV_lesions)

UVLesion_MT<-cbind(UVLesion_MT, ranked_UV_lesions)
colnames(UVLesion_MT)[c(8,9)]<-c("ranked_CPD","ranked_64PP")

# Convert dataframe to a GRanges type and sort according to chromosome order (e.g. chr1, chr2, ect)
UVLesion_MT<-makeGRangesFromDataFrame(UVLesion_MT, keep.extra.columns = T)
UVLesion_MT<-sortSeqlevels(UVLesion_MT)
UVLesion_MT<-sort(UVLesion_MT)
isSorted(UVLesion_MT)

# Measure the difference between ranked 6-4 and CPD
Ranked_diff<-UVLesion_MT
Ranked_diff$Difference<-Ranked_diff$ranked_64PP-Ranked_diff$ranked_CPD

# Change back to data.frame and save
Ranked_diff<-as.data.frame(Ranked_diff)
save(Ranked_diff, file = paste0(out, "UV lesions rank normalized to norm dist 1Mb.RData"))






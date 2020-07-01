library(GenomicRanges)
library(rtracklayer)


## Description: This script measures the mean fold change signal for all IMR90 epigenetic features (E017) binned within chromatin states

out<-"output_directory/"


# Load IMR90 15 core chromatin state dataset
core15<-file.path("file_path/E017_15_coreMarks_dense.bed.bgz")
Chrom15<-import.bed(core15)
mcols(Chrom15)$score<-NULL
mcols(Chrom15)$itemRgb<-NULL
mcols(Chrom15)$thick<-NULL

#

pathToBigWig<-file.path("path_to_bigwig_histone_files//")
HistModFiles<-list.files("path_to_bigwig_histone_files//")

Split<-unlist(strsplit(HistModFiles, "E017-"))
Split<-unlist(strsplit(Split, ".fc.signal.bigwig"))

alldf = as.data.frame(matrix(nrow=16, ncol=length(HistModFiles)))
colnames(alldf) = HistModFiles
chrom_state<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
rownames(alldf)<- c(chrom_state, "genome mean")

for (j in 1:length(HistModFiles)) {
  bigwig<- paste0(pathToBigWig,HistModFiles[j] )
  mark <- import.bw(con=bigwig, as='GRanges')
  mark$score[which(mark$score==0)] <- NA
  mark<-as.data.frame(mark)
  mark = na.omit(mark)
  markmed = mean(na.omit(mark$score))
  mark$score = mark$score/markmed
  alldf[16,j] =  mean(mark$score)
  mark<-makeGRangesFromDataFrame(mark,keep.extra.columns = T)
  
  for (i in 1:length(chrom_state)) {
    state<-as.data.frame(Chrom15)
    state<-state[which(state$name==chrom_state[i]),]
    state<-makeGRangesFromDataFrame(state, keep.extra.columns = T)
    statesig <- subsetByOverlaps(mark,state)
    alldf[i,j] =  mean(statesig$score)
    
  }
  write.csv(alldf, paste0(out, "E017_summary_StateMean_GenomeMean.csv"))
  rm(mark)
  rm(state)
}
alldf


## END ##
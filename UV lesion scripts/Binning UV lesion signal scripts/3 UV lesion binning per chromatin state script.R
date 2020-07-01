library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Description: This is the script used to bin 6-4PP and CPD signal track files (bigwig format) to chromatin state bin sizes (various bin sizes)


# E017_15_coreMarks_dense.bed.bgz file was downloaded from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/


out<-"output_directory/"

# Import core15 bed file
core15<-"path_to_IMR90_chromatin_state_BED_file/E017_15_coreMarks_dense.bed.bgz"
df_core15<-import.bed(core15)
mcols(df_core15)$score<-NULL
mcols(df_core15)$itemRgb<-NULL
mcols(df_core15)$thick<-NULL
df_core15


######
## Import IP64 bw file as RleList
IP64bw<-"path_to_IP64_signal_track_file/IP64100vInput.fc.signal.bw"
dat <- import.bw(con=IP64bw, as='RleList')
dat[dat == 0]<- NA

seqlevels(df_core15)<-seqlevels(dat)
df_core15<- sort(df_core15)

IP64_15<- binnedAverage(bins = df_core15, numvar = dat, varname = "IP64_mean_FC", na.rm = TRUE) # Slow, takes about 1.5-2 hrs
IP64_15<-sortSeqlevels(IP64_15)
IP64_15<-sort(IP64_15)

save(IP64_15, file= paste0(out, "IP64_binned_to_chrom15.RData"))


######
## Import CPD bw file as RleList
CPD_bw<-"path_to_CPD_signal_track_file/CPD100vInput.fc.signal.bw"
dat <- import.bw(con=CPD_bw, as='RleList')
dat[dat == 0]<- NA

seqlevels(df_core15)<-seqlevels(dat)
df_core15<- sort(df_core15)

CPD_15<- binnedAverage(bins = df_core15, numvar = dat, varname = "CPD_mean_FC", na.rm = TRUE) # Slow, takes about 1.5-2 hrs
CPD_15<-sortSeqlevels(CPD_15)
CPD_15<-sort(CPD_15)

save(CPD_15, file= paste0(out, "CPD_binned_to_chrom15.RData"))



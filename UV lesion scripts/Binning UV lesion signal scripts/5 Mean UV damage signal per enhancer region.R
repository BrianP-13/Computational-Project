library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


# Description: This script will measure the average UV lesion signal (fold change) over all annotated enhancers for IMR90 cells.

out<-"output_directory/"

## Load IMR90 enhancer dataset
IMR90_Enh_pos<-import.bed(con = 'file_path/IMR90_enhancer_list_EnhancerAtlas2.0.bed')
IMR90_Enh_pos<-as.data.frame(IMR90_Enh_pos)
IMR90_Enh_pos$start<-IMR90_Enh_pos$start-1
IMR90_Enh_pos<-makeGRangesFromDataFrame(IMR90_Enh_pos,keep.extra.columns = T)
genome(IMR90_Enh_pos)<-"hg19"
IMR90_Enh_pos<-sortSeqlevels(IMR90_Enh_pos)
IMR90_Enh_pos<-sort(IMR90_Enh_pos)
IMR90_Enh_pos

seqlevelsStyle(IMR90_Enh_pos)
isSorted(IMR90_Enh_pos)
seqinfo(IMR90_Enh_pos)
#seqinfo(BSgenome.Hsapiens.UCSC.hg19)



######
## 6-4PP

# Binning 6-4PP FC per gene (pooled signal)
IP64_bw<-'path_to_IP64_signal_track_file/IP64100vInput.fc.signal.bw'
IP64 <- import.bw(con= IP64_bw, as='RleList')
IP64[IP64 == 0] <- NA

seqinfo(IP64)<-seqinfo(IMR90_Enh_pos)
seqlevels(IMR90_Enh_pos)<-names(IP64)

IP64_Enh<- binnedAverage(IMR90_Enh_pos,numvar=IP64,varname='Mean_64PP_FC', na.rm = TRUE)
IP64_Enh<-sortSeqlevels(IP64_Enh)
IP64_Enh<-sort(IP64_Enh)

save(IP64_Enh, file= paste0(out,'IP64_FC_per_Enh.RData'))



######

# Binning CPD FC per gene (Pooled replicates)
CPD_bw<-'path_to_CPD_signal_track_file/CPD100vInput.fc.signal.bw'
CPD <- import.bw(con= CPD_bw, as='RleList')
CPD[CPD == 0] <- NA

seqinfo(CPD)<-seqinfo(IMR90_Enh_pos)
seqlevels(IMR90_Enh_pos)<-names(CPD)

CPD_Enh<- binnedAverage(IMR90_Enh_pos,numvar=CPD,varname='Mean_CPD_FC', na.rm = TRUE)
CPD_Enh<-sortSeqlevels(CPD_Enh)
CPD_Enh<-sort(CPD_Enh)
isSorted(CPD_Enh)

save(CPD_Enh, file= paste0(out,'CPD_FC_per_Enh.RData'))


### END ###
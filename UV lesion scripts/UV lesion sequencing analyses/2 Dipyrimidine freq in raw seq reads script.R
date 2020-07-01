library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicAlignments)


## Description: This script measures the di-pyrimidine frequency from raw sequencing reads of Input and IP-6-4PP for both replicates. Final di-pyrimidine frequencies are averages between both replicates
#

out<-"output_directory/"

# Di-nucleotide function
dinucleotideFrequencyAlong <- function(x,dinucs=c('AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT'),step=1,as.prob=FALSE) {
  w <- width(x)
  if(length(table(w))!=1) stop('all elements in x must have the same width')
  w <- unique(w)
  at.start <- seq(1,(w-step),by=step)
  at.end <- at.start+1
  dat <- lapply(seq_along(at.start),function(i) {
    tmp <- melt(nucleotideFrequencyAt(x,at=c(at.start[[i]],at.end[[i]]),as.prob = as.prob))
    tmp$dinuc <- paste0(tmp$Var1,tmp$Var2)
    tmp[tmp$dinuc %in% dinucs,3:4]
  })
  tmpp <- do.call(rbind,dat)
  toRet <- unlist(lapply(seq_along(dinucs),function(i) {
    dinuc.i <- dinucs[[i]]
    val <- list(tmpp[tmpp$dinuc==dinuc.i,1])
    names(val) <- dinuc.i
    val
  }),recursive=FALSE)
  toRet <- do.call(cbind,toRet)
  return(toRet)
}


## file paths to FASTQ files
ctl1<-file.path('path_to_Input_FASTQ_replicate_1/Input0-1-IL5840-11_S9.fastq.gz')
ctl2<-file.path('path_to_Input_FASTQ_replicate_2/Input0-2-IL5840-12_S10.fastq.gz')
ip64_1<-file.path('path_to_IP64_FASTQ_replicate_1/lane2_AGTG_L002_R1_ALL.fastq.gz')
ip64_2<-file.path('path_to_IP64_FASTQ_replicate_2/lane2_CAAT_L002_R1_ALL.fastq.gz')

# Input
FASTQ1<-readDNAStringSet(filepath = ctl1, format = "fastq",
                         nrec=-1L, skip=0L, seek.first.rec=FALSE,
                         use.names=FALSE, with.qualities=FALSE)
FASTQ2<-readDNAStringSet(filepath = ctl2, format = "fastq",
                         nrec=-1L, skip=0L, seek.first.rec=FALSE,
                         use.names=FALSE, with.qualities=FALSE)
FASTQ3 <- c(FASTQ1,FASTQ2) # Combine list of sequence reads together, then measure the di nucleotide frequency
Pooled_Freq<-dinucleotideFrequencyAlong(FASTQ3, as.prob = T)
Pooled_Freq<-as.data.frame(Pooled_Freq)
write.csv(Pooled_Freq, paste0(out, "64PP_pooled_Input_diPy_FREQ_FASTQ.csv", sep=""))

rm(FASTQ1)
rm(FASTQ2)
rm(FASTQ3)


# IP100
FASTQ1<-readDNAStringSet(filepath = ip64_1, format = "fastq",
                         nrec=-1L, skip=0L, seek.first.rec=FALSE,
                         use.names=FALSE, with.qualities=FALSE)
FASTQ2<-readDNAStringSet(filepath = ip64_2, format = "fastq",
                         nrec=-1L, skip=0L, seek.first.rec=FALSE,
                         use.names=FALSE, with.qualities=FALSE)
FASTQ3 <- c(FASTQ1,FASTQ2) # Combine list of sequence reads together, then measure the di nucleotide frequency
Pooled_Freq<-dinucleotideFrequencyAlong(FASTQ3, as.prob = T)
Pooled_Freq<-as.data.frame(Pooled_Freq)
write.csv(Pooled_Freq, paste0(out, "64PP_pooled_IP100_diPy_FREQ_FASTQ.csv", sep=""))

rm(FASTQ1)
rm(FASTQ2)
rm(FASTQ3)


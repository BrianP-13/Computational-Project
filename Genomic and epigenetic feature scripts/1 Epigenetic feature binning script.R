library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Description: This is the script to bin all epigenetic feature signal track files (bigwig format) for IMR90 (E017)

out<-"output_directory/"


args = commandArgs(trailingOnly = TRUE)
pathvar <- args[1]
filenames <- list.files(path=pathvar)
print(length(filenames))
for(i in 1:length(filenames)) {
  
  bw<-paste(pathvar,filenames[i],sep="")
  print(bw)
  
  dat <- import.bw(con=bw, as='RleList')
  dat[dat == 0] <- NA
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  
  si.gen <- seqinfo(gen)
  
  si <- si.gen[names(dat)]
  
  bins <- tileGenome(si,tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  
  write.csv(as.data.frame(ba),paste(out,filenames[i],"_binned_1Mb.csv",sep=""))
  
}

## END ##
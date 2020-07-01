library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

### Description: This is the script used to bin 6-4PP and CPD signal track files (bigwig format) at 100kb and 1Mb bin size.


out<-"output_directory/"

### Bin size = 100kb

pathvar<-'path_to_directory_with_UV_lesion_bigwig_files/'
filenames <- list.files(path=pathvar)

for(i in 1:length(filenames)) {
  
  bw<-paste(pathvar,filenames[i],sep="")
  print(bw)
  
  dat <- import.bw(con=bw, as='RleList')
  dat[dat == 0] <- NA
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  
  si.gen <- seqinfo(gen)
  
  si <- si.gen[names(dat)]
  
  bins <- tileGenome(si,tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)
  
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  
  write.csv(as.data.frame(ba), paste(out,filenames[i], "binned_100kb.csv", sep=""))
  
}


### Bin size= 1MB

pathvar<-'path_to_directory_with_UV_lesion_bigwig_files/'
filenames <- list.files(path=pathvar)

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
  
  write.csv(as.data.frame(ba), paste(out,filenames[i], "binned_1MB.csv", sep=""))
  
}



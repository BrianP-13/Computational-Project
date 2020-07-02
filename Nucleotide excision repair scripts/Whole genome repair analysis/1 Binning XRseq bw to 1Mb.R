library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: This script is the first step in measuring the excision repair levels for 6-4PP and CPD lesions across the genome (1Mb bin size). Here, processed XRseq bigwig files
# are averaged over 1Mb genomic regions. There are two replicates further separated by plus and minus strand for all repair time points.


out<-"output_directory/"


######
# paths to CPD XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files("path_to_XRseq_files/")
XR_Names<-list()
for (i in 1:length(file_names)) {
  Names<-unlist(strsplit(file_names[i], "*_NHF1"))[2]
  Names<-unlist(strsplit(Names, "_UNIQUE_NORM_fixedStep_25.bw"))[1]
  XR_Names<-append(XR_Names, Names)
}
XR_Names<-unlist(XR_Names)

for(i in 1:length(file_names)) {
  bw<-paste(path,file_names[i],sep="/")
  print(XR_Names[i])
  
  dat <- import.bw(con=bw, as='RleList')
  dat[dat == 0] <- NA
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  si.gen <- seqinfo(gen)
  si <- si.gen[names(dat)]
  bins <- tileGenome(si,tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  ba<- data.frame(ba)
  
  save( ba, file = paste0(out,XR_Names[i], "_binned_1Mb.RData"))
}


######
##

# paths to 64PP XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files("path_to_XRseq_files/")

XR_Names<-list()
for (i in c(1:20)) {
  Names<-unlist(strsplit(file_names[i], "*_NHF1"))[2]
  Names<-unlist(strsplit(Names, "_UNIQUE_NORM_fixedStep_25.bw"))[1]
  XR_Names<-append(XR_Names, Names)
}
XR_Names<-unlist(XR_Names)

for(i in 1:length(file_names)) {
  bw<-paste(path,file_names[i],sep="/")
  print(XR_Names[i])
  
  dat <- import.bw(con=bw, as='RleList')
  dat[dat == 0] <- NA
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  si.gen <- seqinfo(gen)
  si <- si.gen[names(dat)]
  bins <- tileGenome(si,tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  ba<- data.frame(ba)
  
  save(ba, file = paste0(out,XR_Names[i], "_binned_1Mb.RData"))
}



## END ##

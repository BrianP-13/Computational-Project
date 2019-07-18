library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# CPD files

out<-'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/DNA Lesion Data/binned_files/BinnedByAverage/CPD_64PP_binned/CPD/Pooled/'

filepath<-file.path('D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/DNA Lesion Data/GSE94434_DAK_CPD100vInput.fc.signal.bw')
binsizes<-c('1000', '5000', '10000', '50000', '100000', '500000', '1000000')
binsizes<-rev(binsizes)
binsizes[7]

binNames<-c('1kb','5kb','10kb','50kb','100kb','500kb','1Mb')
binNames<-rev(binNames)
binNames[7]

dat <- import.bw(con=filepath, as='RleList')
dat[dat == 0] <- NA

gen <- BSgenome.Hsapiens.UCSC.hg19
si.gen <- seqinfo(gen)
si <- si.gen[names(dat)]

i=7
print(paste0("working on ", binNames[i]))
Sys.time()
num<-as.numeric(binsizes[i])
bins <- tileGenome(si, tilewidth = num, cut.last.tile.in.chrom = TRUE)
ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
write.csv(as.data.frame(ba), paste0(out, "IPCPD100_pooled_", binNames[i], ".csv", sep=""))


for (i in 6:7) {
  print(paste0("working on ", binNames[i]))
  num<-as.numeric(binsizes[i])
  bins <- tileGenome(si, tilewidth = num, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, "IPCPD100_pooled_", binNames[i], ".csv", sep=""))
  
}

# 6-4 PP files
out<-'D:/Sequencing Data/Temp_bigwig/binned/'

filepath<-'D:/Sequencing Data/Temp_bigwig/IP64/'
fileNames<-c("lane2_AGTG_L002_R1_ALL.nodup.tagAlign_x_Input0-1-IL5840-11_S9.nodup_pooled.tagAlign.fc.signal.bw", "lane2_CAAT_L002_R1_ALL.nodup.tagAlign_x_Input0-1-IL5840-11_S9.nodup_pooled.tagAlign.fc.signal.bw", "lane2_AGTG_L002_R1_ALL.nodup_pooled.tagAlign_x_Input0-1-IL5840-11_S9.nodup_pooled.tagAlign.fc.signal.bw")
fileName<-c("IP64_rep1_", "IP64_rep2_", "IP64_pooled_")

for (i in 1:3) {
  bw<-paste0(filepath, fileNames[i])
  print(paste0("working on ", fileNames[i]))
  
  dat <- import.bw(con=bw, as='RleList')
  print("bigwig import complete")
  dat[dat == 0] <- NA
  
  gen <- BSgenome.Hsapiens.UCSC.hg19
  
  si.gen <- seqinfo(gen)
  
  si <- si.gen[names(dat)]
  
  print(paste0("setting bin size at 1e3"))
  
  bins <- tileGenome(si, tilewidth = 1e3, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, fileName[i], "1kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e3"))
  
  bins <- tileGenome(si, tilewidth = 5e3, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, fileName[i], "5kb.csv", sep=""))
  
  print(paste0("setting bin size at 1e4"))
  
  bins <- tileGenome(si, tilewidth = 1e4, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, fileName[i], "10kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e4"))
  
  bins <- tileGenome(si, tilewidth = 5e4, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, fileName[i], "50kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e5"))
  
  bins <- tileGenome(si, tilewidth = 5e5, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out, fileName[i], "500kb.csv", sep=""))
  
  
}


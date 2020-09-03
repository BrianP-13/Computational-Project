library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# Description: This script will measure the mean 6-4PP lesion signal for each replicate across the genome to compare Pearson correlation at various bin sizes


######
## Binning replicate 6-4PP signal track files
out<-c("output_directory/Rep1/","output_directory/Rep2/")
filepath<-'path_to_replicate_IP64_bigwig_files'

fileNames<-c("lane2_AGTG_L002_R1_ALL.nodup.tagAlign_x_Input0-1-IL5840-11_S9.nodup_pooled.tagAlign.fc.signal.bw", "lane2_CAAT_L002_R1_ALL.nodup.tagAlign_x_Input0-1-IL5840-11_S9.nodup_pooled.tagAlign.fc.signal.bw")
fileName<-c("IP64_rep1_", "IP64_rep2_")

for (i in 1:2) {
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
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "1kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e3"))
  
  bins <- tileGenome(si, tilewidth = 5e3, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "5kb.csv", sep=""))
  
  print(paste0("setting bin size at 1e4"))
  
  bins <- tileGenome(si, tilewidth = 1e4, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "10kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e4"))
  
  bins <- tileGenome(si, tilewidth = 5e4, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "50kb.csv", sep=""))
  
  print(paste0("setting bin size at 5e5"))
  
  bins <- tileGenome(si, tilewidth = 5e5, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "500kb.csv", sep=""))
  
  print(paste0("setting bin size at 1e6"))
  
  bins <- tileGenome(si, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  ba <- binnedAverage(bins,numvar=dat,varname='average', na.rm = TRUE)
  write.csv(as.data.frame(ba), paste0(out[i], fileName[i], "1Mb.csv", sep=""))
  
}


######
## Pearson correlation test and line graph

Path1<-"output_directory/Rep1/"
Path2<-"output_directory/Rep2/"

binsize<-c("1kb","5kb","10kb","50kb","100kb","500kb","1Mb")
Corr_table<-data.frame(row.names = "pearson correlation")

for (i in 1:7) {
  filepath1<-list.files(Path1, pattern = binsize[i])
  filepath2<-list.files(Path2, pattern = binsize[i])
  df1<-read.csv(paste0(Path1, filepath1), header = TRUE, sep = ",")
  df2<-read.csv(paste0(Path2, filepath2), header = TRUE, sep = ",")
  
  df3<-cbind(df1, df2$average)
  df3<-subset.data.frame(df3, seqnames != "chrY")
  df3<-subset.data.frame(df3, seqnames != "chrM")
  df3<-na.omit(df3)
  
  cor<-cor(df3[,7], df3[,8])
  Corr_table<-cbind(Corr_table, cor)
  colnames(Corr_table)[i]<-binsize[i]
}
Corr_table

## Line plot
df1<-Corr_table
bins<-c("1kb", "5kb", "10kb", "50kb", "100kb", "500kb", "1Mb")

df1<-t(df1)
df1<-as.data.frame(df1)
df1$binsize<-bins

colnames(df1)[1]<-"Pearson Correlation"
df1$binsize<-factor(x=df1$binsize, levels = c("1kb", "5kb", "10kb", "50kb", "100kb", "500kb", "1Mb"))

# Line plot with labels
IP64_line<-ggplot(df1, aes(x=df1[,2], y=df1[,1], group=1))+geom_line()+geom_point()+ 
  ylim(0.55, 1.0) + xlab("bin size") + ylab("Pearson Correlation")+
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8),
                   panel.grid.minor = element_blank())
IP64_line
ggsave(filename = "64PP_rep_corr_plot.pdf",
       path = out,plot = IP64_line,device = "pdf",width = 2,height = 2,units = "in")

# Line plot - no labels
IP64_line2<-ggplot(df1, aes(x=df1[,2], y=df1[,1], group=1))+geom_line(color="grey")+geom_point()+ 
  ylim(0.55, 1.0) + xlab("bin size") + ylab("Pearson Correlation")+
  theme_bw()+ theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    panel.grid.minor = element_blank())
IP64_line2
ggsave(filename = "64PP_rep_corr_plot.NO.LABS.pdf",
       path = out,plot = IP64_line2,device = "pdf",width = 2,height = 1.5,units = "in")




### END ###

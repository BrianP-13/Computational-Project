library(GenomicRanges)
library(ggplot2)

## Description: Scatter plot analysis comparing genome-wide cumulative nucleotide excision repair and rank normalized UV lesion abundances

## Strategy:
# step 1: Load final 1Mb XRseq data table for each UV lesion type.
# step 2: Calculate the cumulative NER by adding together the NER 
# step 3: Load table with rank normalized UV lesion values. Table with rank normalized values 
# will have some rows missing so you will have to match the overlapping regions from each data set.
# step 4: Conduct Pearson's correlation test to determine the degree of 
# correlation. Make sure to remove chromosomes X,Y, and M from the analysis.
# step 5: Create scatter plot.


out<-"output_directory/"


# Stats table
Stats_table<-as.data.frame(matrix(nrow = 4,ncol = 3))
colnames(Stats_table)<-c("Name","Pearsons_corr","PVAL")
Stats_table[,1]<-c("64PP_FC","64PP_Ranked","CPD_FC","CPD_Ranked")

# Load rank normalized dataset
load(file = "file_path/UV lesions rank normalized to norm dist.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff,keep.extra.columns = T)


######
## 6-4PP

# Load XRseq dataset
load(file = "file_path/IP64_lesion_XRseq_mastertable.RData")
IP64_MT$Total_64NER<-rowSums(x = IP64_MT[,c(7:11)])
IP64_MT<-IP64_MT[,c(1:6,12)]

IP64_MT<-makeGRangesFromDataFrame(IP64_MT,keep.extra.columns = T)
IP64_MT<-sortSeqlevels(IP64_MT)
IP64_MT<-sort(IP64_MT)

# Combine ranked UV lesions signals with total NER data
XR64_MT<-subsetByOverlaps(x = IP64_MT,ranges = Ranked_diff)
mcols(XR64_MT)<-cbind(mcols(XR64_MT), mcols(Ranked_diff)[4])
mcols(XR64_MT)<-mcols(XR64_MT)[c(1,3,2)]


## 64PP fold change vs total 64PP NER

# Process dataset: Remove NAs and log2 transform 64PP FC values
XR64_v_64FC<-as.data.frame(XR64_MT)
XR64_v_64FC<-subset(XR64_v_64FC,seqnames !="chrX")
XR64_v_64FC<-XR64_v_64FC[,c(6,8)]
XR64_v_64FC[,1]<-log2(XR64_v_64FC[,1])
XR64_v_64FC<-na.omit(XR64_v_64FC)

# Correlation test statistic
CorTest<-cor.test(XR64_v_64FC[,1],XR64_v_64FC[,2])

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))
if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}
CorTest

Stats_table[1,2]<-CorTest$estimate
Stats_table[1,3]<-CorTest$p.value

xlab<-"6-4PP Lesions log2(FC)"
ylab<-"Total 6-4PP repair "

## Scatter plot with labels
XR64_v_64FC_plot<-ggplot(XR64_v_64FC, aes(x=XR64_v_64FC[,1], y=XR64_v_64FC[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  labs(x= xlab, y= ylab) +
  annotate(geom="text",label=corr_coef,size=4,fontface="bold.italic", x=quantile(XR64_v_64FC[,1], 0.999), y=quantile(range(min(quantile(XR64_v_64FC[,2],0.005) ), max(XR64_v_64FC[,2])), 0.8),
           hjust=0, vjust=-1) +
  annotate(geom="text",label=pval_coef,size=4,  fontface="bold.italic", x=quantile(XR64_v_64FC[,1], 0.999), y=quantile(range(min(quantile(XR64_v_64FC[,2],0.005) ), max(XR64_v_64FC[,2])), 0.8),
           hjust=0, vjust=0.5)+
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
XR64_v_64FC_plot
ggsave(filename = "IP64PP_FC_vs_sumNER_1Mb.pdf",
       path = out,plot = XR64_v_64FC_plot,device = "pdf",width = 2,height = 2,units = "in")



## Ranked 64PP vs total 64PP NER

# Process dataset: Remove NAs and log2 transform 64PP FC values
XR64_v_ranked64<-as.data.frame(XR64_MT)
XR64_v_ranked64<-subset(XR64_v_ranked64,seqnames !="chrX")
XR64_v_ranked64<-XR64_v_ranked64[,c(7,8)]
XR64_v_ranked64<-na.omit(XR64_v_ranked64)

# Correlation test statistic
CorTest<-cor.test(XR64_v_ranked64[,1],XR64_v_ranked64[,2])

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))
if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

Stats_table[2,2]<-CorTest$estimate
Stats_table[2,3]<-CorTest$p.value

xlab<-"6-4PP Lesions (Rank normalized)"
ylab<-"Total 6-4PP repair "

## Scatter plot with labels
XR64_v_Ranked64_plot<-ggplot(XR64_v_ranked64, aes(x=XR64_v_ranked64[,1], y=XR64_v_ranked64[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  labs(x= xlab, y= ylab) +
  annotate(geom="text",label=corr_coef,size=4,fontface="bold.italic", x=quantile(XR64_v_ranked64[,1], 0.9), y=quantile(range(min(quantile(XR64_v_ranked64[,2],0.005) ), max(XR64_v_ranked64[,2])), 0.8),
           hjust=0, vjust=-1) +
  annotate(geom="text",label=pval_coef,size=4,  fontface="bold.italic", x=quantile(XR64_v_ranked64[,1], 0.9), y=quantile(range(min(quantile(XR64_v_ranked64[,2],0.005) ), max(XR64_v_ranked64[,2])), 0.8),
           hjust=0, vjust=0.5)+
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
XR64_v_Ranked64_plot
ggsave(filename = "IP64PP_Ranked_vs_sumNER_1Mb.pdf",
       path = out,plot = XR64_v_Ranked64_plot,device = "pdf",width = 2,height = 2,units = "in")


######
## CPD damage vs total CPD NER

# Load XRseq dataset
load(file = "file_path/CPD_lesion_XRseq_mastertable.RData")
CPD_MT$Total_cpdNER<-rowSums(x = CPD_MT[,c(7:12)])
CPD_MT<-CPD_MT[,c(1:6,13)]

CPD_MT<-makeGRangesFromDataFrame(CPD_MT,keep.extra.columns = T)
CPD_MT<-sortSeqlevels(CPD_MT)
CPD_MT<-sort(CPD_MT)

# Combine ranked UV lesions signals with total NER data
XRcpd_MT<-subsetByOverlaps(x = CPD_MT,ranges = Ranked_diff)
mcols(XRcpd_MT)<-cbind(mcols(XRcpd_MT), mcols(Ranked_diff)[3])
mcols(XRcpd_MT)<-mcols(XRcpd_MT)[c(1,3,2)]


## CPD fold change vs total CPD NER

# Process dataset: Remove NAs and log2 transform CPD FC values
XRcpd_v_cpdFC<-as.data.frame(XRcpd_MT)
XRcpd_v_cpdFC<-subset(XRcpd_v_cpdFC,seqnames !="chrX")
XRcpd_v_cpdFC<-XRcpd_v_cpdFC[,c(6,8)]
XRcpd_v_cpdFC[,1]<-log2(XRcpd_v_cpdFC[,1])
XRcpd_v_cpdFC<-na.omit(XRcpd_v_cpdFC)

# Correlation test statistic
CorTest<-cor.test(XRcpd_v_cpdFC[,1],XRcpd_v_cpdFC[,2])

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))
if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

Stats_table[3,2]<-CorTest$estimate
Stats_table[3,3]<-CorTest$p.value

xlab<-"CPD Lesions log2(FC)"
ylab<-"Total CPD repair "

## Scatter plot with labels
XRcpd_v_cpdFC_plot<-ggplot(XRcpd_v_cpdFC, aes(x=XRcpd_v_cpdFC[,1], y=XRcpd_v_cpdFC[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  labs(x= xlab, y= ylab) +
  annotate(geom="text",label=corr_coef,size=4,fontface="bold.italic", x=quantile(XRcpd_v_cpdFC[,1], 0.2), y=quantile(range(min(quantile(XRcpd_v_cpdFC[,2],0.5) ), max(XRcpd_v_cpdFC[,2])), 0.8),
           hjust=0, vjust=-1) +
  annotate(geom="text",label=pval_coef,size=4,  fontface="bold.italic", x=quantile(XRcpd_v_cpdFC[,1], 0.2), y=quantile(range(min(quantile(XRcpd_v_cpdFC[,2],0.5) ), max(XRcpd_v_cpdFC[,2])), 0.8),
           hjust=0, vjust=0.5)+
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
XRcpd_v_cpdFC_plot
ggsave(filename = "CPD_FC_vs_sumNER_1Mb.pdf",
       path = out,plot = XRcpd_v_cpdFC_plot,device = "pdf",width = 2,height = 2,units = "in")


## Ranked CPD vs total CPD NER

# Process dataset: Remove NAs and log2 transform 64PP FC values
XRcpd_v_rankedCPD<-as.data.frame(XRcpd_MT)
XRcpd_v_rankedCPD<-subset(XRcpd_v_rankedCPD,seqnames !="chrX")
XRcpd_v_rankedCPD<-XRcpd_v_rankedCPD[,c(7,8)]
XRcpd_v_rankedCPD<-na.omit(XRcpd_v_rankedCPD)

# Correlation test statistic
CorTest<-cor.test(XRcpd_v_rankedCPD[,1],XRcpd_v_rankedCPD[,2])

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))
if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

Stats_table[4,2]<-CorTest$estimate
Stats_table[4,3]<-CorTest$p.value

xlab<-"CPD Lesions (Rank normalized)"
ylab<-"Total CPD repair "

## Scatter plot with labels
XRcpd_v_Rankedcpd_plot<-ggplot(XRcpd_v_rankedCPD, aes(x=XRcpd_v_rankedCPD[,1], y=XRcpd_v_rankedCPD[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  labs(x= xlab, y= ylab) +
  annotate(geom="text",label=corr_coef,size=4,fontface="bold.italic", x=quantile(XRcpd_v_rankedCPD[,1], 0.7), y=quantile(range(min(quantile(XRcpd_v_rankedCPD[,2],0.5) ), max(XRcpd_v_rankedCPD[,2])), 0.8),
           hjust=0, vjust=-1) +
  annotate(geom="text",label=pval_coef,size=4,  fontface="bold.italic", x=quantile(XRcpd_v_rankedCPD[,1], 0.7), y=quantile(range(min(quantile(XRcpd_v_rankedCPD[,2],0.5) ), max(XRcpd_v_rankedCPD[,2])), 0.8),
           hjust=0, vjust=0.5)+
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
XRcpd_v_Rankedcpd_plot
ggsave(filename = "CPD_Ranked_vs_sumNER_1Mb.pdf",
       path = out,plot = XRcpd_v_Rankedcpd_plot,device = "pdf",width = 2,height = 2,units = "in")



write.csv(x = Stats_table,file = paste0(out,"Scatterplot_stats.csv"))



## END ##
library(GenomicRanges)
library(ggplot2)
library(ggrepel)


# Description: This script conducts Principle Component Analysis (PCA) using all epigenetic features and UV lesion signals binned at 1Mb


out<-"output_directory/"


## Import Master Table (Fold Change), does not contain rank normalized UV lesion signals
load(file = "file_path/Master_Table_FC_1Mb.No_RankNorm.RData")
df_MT<-makeGRangesFromDataFrame(df_MT,keep.extra.columns = T)

df_MT<-subset(df_MT, seqnames!="chrY")
df_MT<-subset(df_MT, seqnames!="chrM")
seqlevels(df_MT)<-seqlevelsInUse(df_MT)

df_MT<-as.data.frame(df_MT)

# remove certain columns from dataset (genomic ranges info)
PCA_table<-df_MT
PCA_table<-PCA_table[,6:44]

# samples to exclude: LaminB1 FC (col 38), Repli_seq (not wavelet smoothed, col 33)
PCA_table<-PCA_table[,c(-33,-38)]

## Next, log transform dataset (exclude samples that do not require log2 transformation)
# cols to exclude: 33, 37
PCA_table<-na.omit(PCA_table)
PCA_table[,c(-33,-37)]<-log2(PCA_table[,c(-33,-37)])

# Adjust sample names
colnames(PCA_table)[2]<-"6-4PP"
colnames(PCA_table)[33]<-"Repli Seq"
colnames(PCA_table)[34]<-"LMNA"
colnames(PCA_table)[36]<-"DNA methylation"
colnames(PCA_table)[37]<-"LaminB1 (DamID)"

# Conduct PCA analysis
PCA<-prcomp(PCA_table, center = T, scale. = T,rank. = 2, retx = T)
Summary<-summary(PCA)
Summary$importance
SumCor<-as.data.frame(Summary$importance)
#write.csv(SumCor,file = paste0(out, "Summary stats of PCA 1Mb.csv"))

PCAdata<-PCA$x
PCAdata<-as.data.frame(PCAdata)

PCcor<-as.data.frame(matrix(ncol = 2,nrow = 37))  # rows can change depending on input values
colnames(PCcor)<-c("PC1","PC2")
rownames(PCcor)<-colnames(PCA_table)


# Test for loop to calculate the correlation values of each PC to the original dataset
for (i in 1:37) {
  df<-PCA_table[,i]
  PCcor[i,1]<-cor(df, PCAdata[,1])
  PCcor[i,2]<-cor(df, PCAdata[,2])
}


# plot with labels
PCA_plot<-ggplot(data = PCcor, mapping = aes(x=PCcor[,1],y=PCcor[,2], label=rownames(PCcor)))+ 
  geom_point(size=0.25)+
  labs(x=paste0("PC1 Corr. (var. explained = ",signif(SumCor[2,1], digits = 2),")", sep=""),
       y=paste0("PC2 Corr. (var. explained = ",signif(SumCor[2,2], digits = 2),")", sep=""))+
  coord_cartesian(xlim = c(-1.75,1.75), ylim = c(-1,1))+
  geom_vline(xintercept = 0,color="grey",size=0.25)+geom_hline(yintercept = 0,color="grey",size=0.25)+
  geom_text_repel(size=4, aes(label=rownames(PCcor)))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 5),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 5))
PCA_plot
ggsave(filename = "PCA_UV_lesions_and_chromatin_marks.pdf",
       path = out,plot = PCA_plot,device = "pdf",width = 2.5,height = 2.5,units = "in")



## END ##
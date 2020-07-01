library(GenomicRanges)
library(ggplot2)


# Description: This script measures the correlation between 6-4PP lesions and epigenetic features at 1Mb genome size


out<-"output_directory/"


## load master table
load(file = "file_path/Master_Table_FC_1Mb.RData")
df_MT<-subset(df_MT,seqnames != "chrY")
df_MT<-subset(df_MT,seqnames != "chrM")

df_MT2<-df_MT[,c(1:5,7,11,12,28,25,35,24,40,45)]

# Create table to store correlation values and p-values
ColNames<-c('Feature_1','Feature_2', 'Pearsons_r_value','p_value')
Histone_marks<-c("DNase","H2AZ","H3K4me2","H3K36me3","H3K9me3","H3K27me3","H4K20me3","DNAme")

Stats_stable<-as.data.frame(matrix(nrow = 8, ncol = 4))
colnames(Stats_stable)<-ColNames
Stats_stable[,1]<-"6-4PP"
Stats_stable[,2]<-Histone_marks

for (i in 1:8) {
  df<-df_MT2[,c(6,6+i)]
  df<-log2(df)
  df<-na.omit(df)
  CorTest<-cor.test(df[,1], df[,2])
  Stats_stable[i,3]<-CorTest$estimate
  Stats_stable[i,4]<-CorTest$p.value
}
write.csv(Stats_stable,file = paste0(out, "IP64_vs_hist_marks_corr_vals.csv"))


Xlab<-"6-4PP Lesions log2(FC)"

######
## 6-4PP vs DNase
df<-df_MT2[,c(6,7)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"DNase log2(FC)"
IP64_DNase<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$DNase)$stats[1],boxplot.stats(x = df$DNase)$stats[5]),
                     minor_breaks = NULL,breaks = c(0,1,2), expand = expand_scale(add = 0.2))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
        panel.border = element_rect(fill = NA,size = 0.25),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_DNase
ggsave(filename = "64PP_vs_DNase_1Mb.pdf",
       path = out,plot = IP64_DNase,device = "pdf",width = 2,height = 2,units = "in")



######
## 6-4PP vs H2A.Z
df<-df_MT2[,c(6,8)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H2A.Z log2(FC)"
IP64_H2AZ<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H2A.Z)$stats[1],boxplot.stats(x = df$H2A.Z)$stats[5]),
                     minor_breaks = NULL, expand = expand_scale(add = 0.2))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H2AZ
ggsave(filename = "64PP_vs_H2AZ_1Mb.pdf",
       path = out,plot = IP64_H2AZ,device = "pdf",width = 2,height = 2,units = "in")


######
## 6-4PP vs H3K4me2
df<-df_MT2[,c(6,9)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H3K4me2 log2(FC)"
IP64_H3K4me2<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H3K4me2)$stats[1],boxplot.stats(x = df$H3K4me2)$stats[5]),
                     minor_breaks = NULL,breaks = c(-1,0,1), expand = expand_scale(add = 0.25))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K4me2
ggsave(filename = "64PP_vs_H3K4me2_1Mb.pdf",
       path = out,plot = IP64_H3K4me2,device = "pdf",width = 2,height = 2,units = "in")


######
## 6-4PP vs H3K36me3
df<-df_MT2[,c(6,10)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H3K36me3 log2(FC)"
IP64_H3K36me3<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H3K36me3)$stats[1],boxplot.stats(x = df$H3K36me3)$stats[5]),
                     minor_breaks = NULL, expand = expand_scale(add = 0.25))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K36me3
ggsave(filename = "64PP_vs_H3K36me3_1Mb.pdf",
       path = out,plot = IP64_H3K36me3,device = "pdf",width = 2,height = 2,units = "in")



######
## 6-4PP vs H3K9me3
df<-df_MT2[,c(6,11)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H3K9me3 log2(FC)"
IP64_H3K9me3<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H3K9me3)$stats[1],boxplot.stats(x = df$H3K9me3)$stats[5]),
                     minor_breaks = NULL,breaks = c(0,1,2),expand = expand_scale(add = 0.2))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K9me3
ggsave(filename = "64PP_vs_H3K9me3_1Mb.pdf",
       path = out,plot = IP64_H3K9me3,device = "pdf",width = 2,height = 2,units = "in")



######
## 6-4PP vs H3K27me3
df<-df_MT2[,c(6,12)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H3K27me3 log2(FC)"
IP64_H3K27me3<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H3K27me3)$stats[1],boxplot.stats(x = df$H3K27me3)$stats[5]),
                     minor_breaks = NULL,breaks = c(-1,0,1),expand = expand_scale(add = 0.2))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K27me3
ggsave(filename = "64PP_vs_H3K27me3_1Mb.pdf",
       path = out,plot = IP64_H3K27me3,device = "pdf",width = 2,height = 2,units = "in")



######
## 6-4PP vs H4K20me3
df<-df_MT2[,c(6,13)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"H4K20me3 log2(FC)"
IP64_H4K20me3<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$H4K20me3)$stats[1],boxplot.stats(x = df$H4K20me3)$stats[5]),
                     minor_breaks = NULL,expand = expand_scale(add = 0.1))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H4K20me3
ggsave(filename = "64PP_vs_H4K20me3_1Mb.pdf",
       path = out,plot = IP64_H4K20me3,device = "pdf",width = 2,height = 2,units = "in")



######
## 6-4PP vs DNAme
df<-df_MT2[,c(6,14)]
df<-log2(df)
df<-na.omit(df)

Ylab<-"DNAme log2(FC)"
IP64_DNAme<-ggplot(df, aes(x=df[,1], y=df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = df$IP64)$stats[1],boxplot.stats(x = df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = df$WGBS_fractional_methylation)$stats[1],boxplot.stats(x = df$WGBS_fractional_methylation)$stats[5]),
                     minor_breaks = NULL,expand = expand_scale(add = 0.2))+
  labs(x= Xlab, y= Ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_DNAme
ggsave(filename = "64PP_vs_DNAme_1Mb.pdf",
       path = out,plot = IP64_DNAme,device = "pdf",width = 2,height = 2,units = "in")


## END ##
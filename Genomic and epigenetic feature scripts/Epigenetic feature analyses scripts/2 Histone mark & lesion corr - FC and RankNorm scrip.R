library(GenomicRanges)
library(ggplot2)

# Description: This script measures the Pearson correlation between selected histone marks and UV lesion abundance using fold change (IP/Input) vs rank normalized signal


out<-"output_directory/"


## load master table
load(file = "file_path/Master_Table_FC_1Mb.RData")
df_MT<-subset(df_MT,seqnames != "chrY")
df_MT<-subset(df_MT,seqnames != "chrM")

# Create table to store correlation values and p-values
ColNames<-c('Feature_1','Feature_2', 'Pearsons_r_value','p_value')

Stats_stable<-as.data.frame(matrix(nrow = 4, ncol = 4))
colnames(Stats_stable)<-ColNames
Stats_stable[c(1,3),1]<-"FC_64PP"
Stats_stable[c(2,4),1]<-"Ranked_64PP"
Stats_stable[1,2]<-"FC_CPD"
Stats_stable[2,2]<-"Ranked_CPD"
Stats_stable[3:4,2]<-"FC_H3K9me3"

######
## 6-4PP vs CPD correlation analysis
CPD_64PP_df<-df_MT[,c(1:9)]

## 6-4PP vs CPD (Fold Change)
FC_df<-CPD_64PP_df[,6:7]
FC_df<-log2(FC_df)
FC_df<-na.omit(FC_df)

CorTest<-cor.test(FC_df[,1], FC_df[,2])

Stats_stable[1,3]<-CorTest$estimate
Stats_stable[1,4]<-CorTest$p.value

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))
if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

# scatterplot with labels
xlab<-"CPD Lesions log2(FC)"
ylab<-"6-4PP Lesions log2(FC)"

CPD_64PP_FC_plot<-ggplot(FC_df, aes(x=FC_df[,1], y=FC_df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = FC_df$CPD)$stats[1],boxplot.stats(x = FC_df$CPD)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.8,-0.4,0), expand = expand_scale(add = 0.1))+
  scale_y_continuous(limits = c(boxplot.stats(x = FC_df$IP64)$stats[1],boxplot.stats(x = FC_df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0), expand = expand_scale(add = 0.09))+
  labs(x= xlab, y= ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
        panel.border = element_rect(fill = NA,size = 0.25),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
CPD_64PP_FC_plot
ggsave(filename = "CPD_vs_64PP_FC_1Mb.pdf",
       path = out,plot = CPD_64PP_FC_plot,device = "pdf",width = 2,height = 2,units = "in")


## 6-4PP vs CPD (Rank Normalized)
RankNorm_df<-CPD_64PP_df[,8:9]
RankNorm_df<-na.omit(RankNorm_df)

CorTest<-cor.test(RankNorm_df[,1], RankNorm_df[,2])
Stats_stable[2,3]<-CorTest$estimate
Stats_stable[2,4]<-CorTest$p.value

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))

if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

# scatterplot with labels
xlab<-"CPD Lesions (RankNorm)"
ylab<-"6-4PP Lesions (RankNorm)"

CPD_64PP_Ranked_plot<-ggplot(RankNorm_df, aes(x=RankNorm_df[,1], y=RankNorm_df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = RankNorm_df$ranked_CPD)$stats[1],boxplot.stats(x = RankNorm_df$ranked_CPD)$stats[5]),
                     minor_breaks = NULL,breaks = c(-2,0,2), expand = expand_scale(add = 0.25))+
  scale_y_continuous(limits = c(boxplot.stats(x = RankNorm_df$ranked_64PP)$stats[1],boxplot.stats(x = RankNorm_df$ranked_64PP)$stats[5]),
                     minor_breaks = NULL,breaks = c(-2,0,2), expand = expand_scale(add = 0.25))+
  labs(x= xlab, y= ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
CPD_64PP_Ranked_plot
ggsave(filename = "CPD_vs_64PP_RankNorm_1Mb.pdf",
       path = out,plot = CPD_64PP_Ranked_plot,device = "pdf",width = 2,height = 2,units = "in")


######
## 6-4PP vs H3K9me3
IP64_H3K9_df<-df_MT[,c(1:5,7,9,35)]

# 6-4PP vs H3K9me3 (Fold Change)
FC_df<-IP64_H3K9_df[,c(6,8)]
FC_df<-log2(FC_df)
FC_df<-na.omit(FC_df)

CorTest<-cor.test(FC_df[,1], FC_df[,2])

Stats_stable[3,3]<-CorTest$estimate
Stats_stable[3,4]<-CorTest$p.value

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))

if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

# scatterplot with labels
xlab<-"6-4PP Lesions log2(FC)"
ylab<-"H3K9me3 log2(FC)"

IP64_H3K9_FC_plot<-ggplot(FC_df, aes(x=FC_df[,1], y=FC_df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = FC_df$IP64)$stats[1],boxplot.stats(x = FC_df$IP64)$stats[5]),
                     minor_breaks = NULL,breaks = c(-0.4,-0.2,0),expand = expand_scale(add = 0.05))+
  scale_y_continuous(limits = c(boxplot.stats(x = FC_df$H3K9me3)$stats[1],boxplot.stats(x = FC_df$H3K9me3)$stats[5]),
                     minor_breaks = NULL,expand = expand_scale(add = 0.2))+
  labs(x= xlab, y= ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K9_FC_plot
ggsave(filename = "IP64_H3K9me3_FC_1Mb.pdf",
       path = out,plot = IP64_H3K9_FC_plot,device = "pdf",width = 2,height = 2,units = "in")


# 6-4PP vs H3K9me3 (Rank Normalized)
RankNorm_df<-IP64_H3K9_df[,c(7:8)]
RankNorm_df[,2]<-log2(RankNorm_df[,2])
RankNorm_df<-na.omit(RankNorm_df)

CorTest<-cor.test(RankNorm_df[,1], RankNorm_df[,2])
Stats_stable[4,3]<-CorTest$estimate
Stats_stable[4,4]<-CorTest$p.value

corr_coef<- paste0("r = ", print(signif(CorTest$estimate, digits=2)))

if (CorTest$p.value < 2.2e-16) {
  pval_coef<- paste0("p < 2.2e-16")
} else {
  pval_coef<- paste0("p = ", print(signif(CorTest$p.value, digits=3)))
}

# scatterplot with labels
xlab<-"6-4PP Lesions (RankNorm)"
ylab<-"H3K9me3 log2(FC)"

IP64_H3K9_Ranked_plot<-ggplot(RankNorm_df, aes(x=RankNorm_df[,1], y=RankNorm_df[,2])) + geom_point(alpha=0.25, size=0.25, shape=19)+ 
  scale_x_continuous(limits = c(boxplot.stats(x = RankNorm_df$ranked_64PP)$stats[1],boxplot.stats(x = RankNorm_df$ranked_64PP)$stats[5]),
                     minor_breaks = NULL,breaks = c(-2,0,2),expand = expand_scale(add = 0.25))+
  scale_y_continuous(limits = c(boxplot.stats(x = RankNorm_df$H3K9me3)$stats[1],boxplot.stats(x = RankNorm_df$H3K9me3)$stats[5]),
                     minor_breaks = NULL,expand = expand_scale(add = 0.25))+
  labs(x= xlab, y= ylab) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8)) +
  geom_smooth(colour= "red", method=lm,size=0.25)
IP64_H3K9_Ranked_plot
ggsave(filename = "64PP_RankNorm_vs_H3K9me3FC_1Mb.pdf",
       path = out,plot = IP64_H3K9_Ranked_plot,device = "pdf",width = 2,height = 2,units = "in")



write.csv(x = Stats_stable,file = paste0(out,"Stats_table_UV_FC_vs_Ranked.csv"))


## END ##
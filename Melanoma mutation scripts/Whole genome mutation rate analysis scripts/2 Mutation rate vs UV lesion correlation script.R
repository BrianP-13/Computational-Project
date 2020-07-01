library(GenomicRanges)
library(IDPmisc)
library(ggplot2)


## Description: This script measure the correlation between mutation rates over 1Mb regions and UV lesion abundance over the whole genome. 
# Correlation analysis of mutation rates will only include rates calculated from the plus strand of the reference genome.


out<-"output_directory/"


######
## Data formating

# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist.RData")
UV_dmg<-makeGRangesFromDataFrame(Ranked_diff,keep.extra.columns = T)

# Import hg19_1Mb_DNA_stats file
load(file = "file_path/hg19_1Mb_DNA_stats.RData")
hg19_nt_stats<-makeGRangesFromDataFrame(hg19_nt_stats,keep.extra.columns = T)

# Import mutation rate file
load(file = "file_path/All_SSM_rates_1Mb_hg19.RData")
mcols(SBS_1Mb)[1:12]<-NULL

#
## Subset data to match genomic ranges in UV_dmg data table
hg19_nt_stats<-subsetByOverlaps(x = hg19_nt_stats,ranges = UV_dmg)
SBS_1Mb<-subsetByOverlaps(x = SBS_1Mb,ranges = UV_dmg)

#
## Normalize UV lesion data by TT feq
# subset standard TT frequency
TT_freq<-as.data.frame(mcols(hg19_nt_stats)[16])
TT_freq<-TT_freq[,1]/median(TT_freq[,1])

# Create data.frame for each UV lesion type
IP64<-as.data.frame(UV_dmg)
IP64<-IP64[,c(1:5,7)]
IP64$Norm_64PP<-IP64$IP64*TT_freq
IP64$Norm_64PP<-log2(IP64$Norm_64PP)
IP64<-as.data.frame(cbind(IP64,mcols(SBS_1Mb)[1:12]))

CPD<-as.data.frame(UV_dmg)
CPD<-CPD[,c(1:6)]
CPD$Norm_CPD<-CPD$CPD*TT_freq
CPD$Norm_CPD<-log2(CPD$Norm_CPD)
CPD<-as.data.frame(cbind(CPD, mcols(SBS_1Mb)[1:12]))

# Create table to store correlation values and p-values
ColNames<-c('SBS_type', 'UV_lesion_type', 'Spearmans_rho_value','Spearman_p_value','Pearsons_r_value','Pearson_p_value')

Mut_cor_df<-as.data.frame(matrix(nrow = 2, ncol = 6))
Mut_cor_df[,1]<-"C_to_T"
colnames(Mut_cor_df)<-ColNames
Mut_cor_df[1,2]<-"64PP"
Mut_cor_df[2,2]<-"CPD"


#
## 64PP correlation
df_64<-IP64[,c(7,11)]

# Spearmans R correlation
S_CorTest<-cor.test(df_64[,1], df_64[,2], method = "spearman", exact = F)
Rho_coef<- paste0("rho  == ", as.numeric(print(signif(S_CorTest$estimate, digits=2))))
if (S_CorTest$p.value < 2.2e-16) {
  Rho_pval<- paste0("p < 2.2e-16")
} else {
  Rho_pval<- paste0("p = ", print(signif(S_CorTest$p.value, digits=3)))
}

# Pearsons R correlation
P_CorTest<-cor.test(df_64[,1], df_64[,2])
R_coef<- paste0("r = ", print(signif(P_CorTest$estimate, digits=2)))
if (P_CorTest$p.value < 2.2e-16) {
  R_pval<- paste0("p < 2.2e-16")
} else {
  R_pval<- paste0("p = ", print(signif(P_CorTest$p.value, digits=3)))
}

Mut_cor_df[1,3]<- S_CorTest$estimate # Spearman's rho-value
Mut_cor_df[1,4]<- S_CorTest$p.value # p-value for Spearman's rank corr test
Mut_cor_df[1,5]<- P_CorTest$estimate # Pearson's R coefficient (r-value)
Mut_cor_df[1,6]<- P_CorTest$p.value # p-value for Pearson's corr test


# Scatter plot with labels
MAX<-max(df_64[,2])*1.15

XLab<-"6-4PP signal log2(FC)"
YLab<- "C>T mutation rate"

Plot_64<-ggplot(df_64, aes(x=df_64[,1], y=df_64[,2])) + geom_point(alpha=0.25, size=0.25, shape=19) +
  coord_cartesian( ylim = c(0,MAX )) + 
  labs(x= XLab, y= YLab) +
  theme_bw() +
  theme(line = element_line(size = 0.15),
        panel.border = element_rect(fill = NA,size = 0.25),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  annotate(geom="text",label= Rho_coef ,size=3,fontface="bold", hjust=0, vjust=-1, parse=TRUE,
           x=quantile(df_64[,1],0.005),
           y=quantile(  range( min(df_64[,2]), max(df_64[,2]) ), 0.85))Plot_64

ggsave(filename = "IP64_CtoT_mut_rate_1Mb.pdf",
       path = out,plot = Plot_64,device = "pdf",width = 2,height = 2,units = "in")


#
## CPD correlation
df_CPD<-CPD[,c(7,11)]

# Spearmans R correlation
S_CorTest<-cor.test(df_CPD[,1], df_CPD[,2], method = "spearman", exact = F)
Rho_coef<- paste0("rho  == ", as.numeric(print(signif(S_CorTest$estimate, digits=2))))
if (S_CorTest$p.value < 2.2e-16) {
  Rho_pval<- paste0("p < 2.2e-16")
} else {
  Rho_pval<- paste0("p = ", print(signif(S_CorTest$p.value, digits=3)))
}

# Pearsons R correlation
P_CorTest<-cor.test(df_CPD[,1], df_CPD[,2])
R_coef<- paste0("r = ", print(signif(P_CorTest$estimate, digits=2)))
if (P_CorTest$p.value < 2.2e-16) {
  R_pval<- paste0("p < 2.2e-16")
} else {
  R_pval<- paste0("p = ", print(signif(P_CorTest$p.value, digits=3)))
}

Mut_cor_df[2,3]<- S_CorTest$estimate # Spearman's rho-value
Mut_cor_df[2,4]<- S_CorTest$p.value # p-value for Spearman's rank corr test
Mut_cor_df[2,5]<- P_CorTest$estimate # Pearson's R coefficient (r-value)
Mut_cor_df[2,6]<- P_CorTest$p.value # p-value for Pearson's corr test


# Scatter plot with labels
MAX<-max(df_CPD[,2])*1.15

XLab<-"CPD signal log2(FC)"
YLab<- "C>T mutation rate"

Plot_CPD<-ggplot(df_CPD, aes(x=df_CPD[,1], y=df_CPD[,2])) + geom_point(alpha=0.25, size=0.25, shape=19) +
  coord_cartesian( ylim = c(0,MAX )) + 
  labs(x= XLab, y= YLab) +
  theme_bw() +
  theme(line = element_line(size = 0.15),
        panel.border = element_rect(fill = NA,size = 0.25),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  annotate(geom="text",label= Rho_coef ,size=3,fontface="bold", hjust=0, vjust=-1, parse=TRUE,
           x=quantile(df_CPD[,1],0.005),
           y=quantile(  range( min(df_CPD[,2]), max(df_CPD[,2]) ), 0.85))
Plot_CPD

ggsave(filename = "CPD_CtoT_mut_rate_1Mb.pdf",
       path = out,plot = Plot_CPD,device = "pdf",width = 2,height = 2,units = "in")



## END ##
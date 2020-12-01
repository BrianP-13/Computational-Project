library(GenomicRanges)
library(ggplot2)
library(IDPmisc)
library(reshape2)

# Description: Analysis of mutation rates for enhancers with and without excision repair in top UV susceptible regions


out<-"output_directory/"
yLab<-"Mutation rate \n Log10(C>T rate)"

# Overlap function:
Gene_Ovrlp_fun<-function( Query, Subject,Overlap,Region_name){
  df<-Query
  hits<-findOverlaps(query = Query,subject =  Subject,type = "any",ignore.strand=T,select = "all")
  gr.over<- pintersect(Query[queryHits(hits)],Subject[subjectHits(hits)])
  gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
    sum(width(x)))
  df$overlap<-0
  df$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
  df<-data.frame(df)
  df$percent_overlap<-(df$overlap/df$width)*100
  df<-df[which(df$percent_overlap >= Overlap),]
  df<-makeGRangesFromDataFrame(df,keep.extra.columns = T)
  df$Region<-Region_name
  return(df)
}


# Import rank normalized UV dmg dataset and calculate top and bottom 10% regions
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff,keep.extra.columns = T)

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


# Load enhancer dataset 
load(file = 'file_path/Enh_Total_XR_and_mutR_dataset.RData')
Enh_XR_mutR<-subset(Enh_XR_mutR,seqnames!="chrX")
Enh_XR_mutR$Log10_mutR<-log10(Enh_XR_mutR$Total_C_to_T_rate)

## Subset enhancer regions to include only 64PP repair
Enh_XR64<-Enh_XR_mutR[,c(1:7,9)]
Enh_XR64$Log10_mutR<-log10(Enh_XR64$Total_C_to_T_rate)

## Subset enhancer regions to include only CPD repair
Enh_XRcpd<-Enh_XR_mutR[,c(1:6,8,9)]
Enh_XRcpd$Log10_mutR<-log10(Enh_XRcpd$Total_C_to_T_rate)



# 6-4PP repair analysis - top/bottom 6-4PP regions ------------------------

# Table for stats testing
Enh_XR64_stats<-as.data.frame( matrix(nrow = 4,ncol = 3))
colnames(Enh_XR64_stats)<-c("Group 1","Group 2","Wilcoxon_RST_pval")

# Label groups for stats test
Enh_XR64_stats[c(1,3),1]<-"Top_10%_64PP-no_64PP_repair"

Enh_XR64_stats[1,2]<-"Top_10%_64PP-64PP_repair"
Enh_XR64_stats[4,1]<-"Top_10%_64PP-64PP_repair"

Enh_XR64_stats[2,1]<-"Bottom_10%_64PP-no_64PP_repair"
Enh_XR64_stats[3,2]<-"Bottom_10%_64PP-no_64PP_repair"

Enh_XR64_stats[2,2]<-"Bottom_10%_64PP-64PP_repair"
Enh_XR64_stats[4,2]<-"Bottom_10%_64PP-64PP_repair"


### Overlap of Enhancers in top and bottom 64PP regions

## Top 10%
XR64_T10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_T10,Overlap = 50,Region_name = "Top 10% 64PP")
XR64_T10_Enh$overlap<-NULL
XR64_T10_Enh$percent_overlap<-NULL
XR64_T10_Enh$Region<-NULL
XR64_T10_Enh<-as.data.frame(XR64_T10_Enh)

# Enhancers in top 10% 64PP regions - no 64PP repair
XR64_T10_Enh_no_repair<-XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair==0),]
XR64_T10_Enh_no_repair$Group<-"Top 10% 64PP - no repair"

# Enhancers in top 10% 64PP regions - with 64PP repair
XR64_T10_Enh_w_repair<-XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair!=0),]
XR64_T10_Enh_w_repair$Group<-"Top 10% 64PP - 64PP repair"

## Bottom 10%
XR64_B10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_B10,Overlap = 50,Region_name = "Bottom 10% 64PP")
XR64_B10_Enh$overlap<-NULL
XR64_B10_Enh$percent_overlap<-NULL
XR64_B10_Enh$Region<-NULL
XR64_B10_Enh<-as.data.frame(XR64_B10_Enh)

# Enhancers in bottom 10% 64PP regions - no 64PP repair
XR64_B10_Enh_no_repair<-XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair==0),]
XR64_B10_Enh_no_repair$Group<-"Bottom 10% 64PP - no repair"

# Enhancers in bottom 10% 64PP regions - with 64PP repair
XR64_B10_Enh_w_repair<-XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair!=0),]
XR64_B10_Enh_w_repair$Group<-"Bottom 10% 64PP - 64PP repair"

## Master table with ALL conditions
ALL_XR64_enh<-rbind(XR64_T10_Enh_no_repair,XR64_T10_Enh_w_repair,XR64_B10_Enh_no_repair,XR64_B10_Enh_w_repair)

# Conduct Wilcoxon rank sum test
Enh_XR64_stats[1,3]<-wilcox.test(x = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Top 10% 64PP - no repair"),9],
                                 y = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Top 10% 64PP - 64PP repair"),9], 
                                 alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XR64_stats[2,3]<-wilcox.test(x = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Bottom 10% 64PP - no repair"),9],
                                 y = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Bottom 10% 64PP - 64PP repair"),9], 
                                 alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XR64_stats[3,3]<-wilcox.test(x = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Top 10% 64PP - no repair"),9],
                                 y = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Bottom 10% 64PP - no repair"),9], 
                                 alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XR64_stats[4,3]<-wilcox.test(x = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Top 10% 64PP - 64PP repair"),9],
                                 y = ALL_XR64_enh[which(ALL_XR64_enh$Group=="Bottom 10% 64PP - 64PP repair"),9], 
                                 alternative = "two.sided", paired = F, conf.int = T)$p.value

# Correct for multiple hypothesis testing 
Enh_XR64_stats$Wilcoxon_RST_pval<-as.numeric(Enh_XR64_stats$Wilcoxon_RST_pval)

Enh_XR64_stats$p.adjust.BH<-p.adjust(p = Enh_XR64_stats$Wilcoxon_RST_pval, method = "BH")
write.csv(Enh_XR64_stats,file = paste0(out, "Wilcoxon RST - Enh mutR in top 10 64PP.csv"))

# Add factor levels to labels
ALL_XR64_enh$Group<-factor(ALL_XR64_enh$Group, levels = c("Top 10% 64PP - no repair","Top 10% 64PP - 64PP repair",
                                                          "Bottom 10% 64PP - no repair","Bottom 10% 64PP - 64PP repair" ),
                           labels = c("Top 10% 64PP -\nno repair","Top 10% 64PP -\n64PP repair",
                                      "Bottom 10% 64PP -\nno repair","Bottom 10% 64PP -\n64PP repair" ))

# Plot
Enh_XR64_All_plot<-ggplot(data = ALL_XR64_enh, aes(x = Group,y = Log10_mutR, fill=Group))+theme_bw()+
  geom_boxplot(show.legend = F)+ylab(yLab)+xlab(NULL)+coord_cartesian(ylim = c(0,-3.25))+
  theme(text = element_text(size = 9),axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(filename = "T10_64PP_Enh_MutRate_w_&_wo_XR64_repair.jpeg",
       path = out,plot = Enh_XR64_All_plot,device = "jpeg",width = 3,height = 4,units = "in")


# CPD repair analysis - top/bottom CPD regions ------------------------

# Table for stats testing
Enh_XRcpd_stats<-as.data.frame( matrix(nrow = 4,ncol = 3))
colnames(Enh_XRcpd_stats)<-c("Group 1","Group 2","Wilcoxon_RST_pval")

# Label groups for stats test
Enh_XRcpd_stats[c(1,3),1]<-"Top_10%_CPD-no_CPD_repair"

Enh_XRcpd_stats[1,2]<-"Top_10%_CPD-CPD_repair"
Enh_XRcpd_stats[4,1]<-"Top_10%_CPD-CPD_repair"

Enh_XRcpd_stats[2,1]<-"Bottom_10%_CPD-no_CPD_repair"
Enh_XRcpd_stats[3,2]<-"Bottom_10%_CPD-no_CPD_repair"

Enh_XRcpd_stats[2,2]<-"Bottom_10%_CPD-CPD_repair"
Enh_XRcpd_stats[4,2]<-"Bottom_10%_CPD-CPD_repair"


### Overlap of Enhancers in top and bottom CPD regions

## Top 10%
XRcpd_T10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = CPD_T10,Overlap = 50,Region_name = "Top 10% CPD")
XRcpd_T10_Enh$overlap<-NULL
XRcpd_T10_Enh$percent_overlap<-NULL
XRcpd_T10_Enh$Region<-NULL
XRcpd_T10_Enh<-as.data.frame(XRcpd_T10_Enh)

# Enhancers in top 10% CPD regions - no CPD repair
XRcpd_T10_Enh_no_repair<-XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair==0),]
XRcpd_T10_Enh_no_repair$Group<-"Top 10% CPD - no CPD repair"

# Enhancers in top 10% CPD regions - with CPD repair
XRcpd_T10_Enh_w_repair<-XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair!=0),]
XRcpd_T10_Enh_w_repair$Group<-"Top 10% CPD - CPD repair"

## Bottom 10%
XRcpd_B10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = CPD_B10,Overlap = 50,Region_name = "Bottom 10% CPD")
XRcpd_B10_Enh$overlap<-NULL
XRcpd_B10_Enh$percent_overlap<-NULL
XRcpd_B10_Enh$Region<-NULL
XRcpd_B10_Enh<-as.data.frame(XRcpd_B10_Enh)

# Enhancers in bottom 10% CPD regions - no CPD repair
XRcpd_B10_Enh_no_repair<-XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair==0),]
XRcpd_B10_Enh_no_repair$Group<-"Bottom 10% CPD - no CPD repair"

# Enhancers in bottom 10% CPD regions - with CPD repair
XRcpd_B10_Enh_w_repair<-XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair!=0),]
XRcpd_B10_Enh_w_repair$Group<-"Bottom 10% CPD - CPD repair"

## Master table with ALL conditions
ALL_XRcpd_enh<-rbind(XRcpd_T10_Enh_no_repair,XRcpd_T10_Enh_w_repair,XRcpd_B10_Enh_no_repair,XRcpd_B10_Enh_w_repair)

# Conduct Wilcoxon rank sum test
Enh_XRcpd_stats[1,3]<-wilcox.test(x = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Top 10% CPD - no CPD repair"),9],
                                  y = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Top 10% CPD - CPD repair"),9], 
                                  alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XRcpd_stats[2,3]<-wilcox.test(x = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Bottom 10% CPD - no CPD repair"),9],
                                  y = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Bottom 10% CPD - CPD repair"),9], 
                                  alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XRcpd_stats[3,3]<-wilcox.test(x = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Top 10% CPD - no CPD repair"),9],
                                  y = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Bottom 10% CPD - no CPD repair"),9], 
                                  alternative = "two.sided", paired = F, conf.int = T)$p.value
Enh_XRcpd_stats[4,3]<-wilcox.test(x = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Top 10% CPD - CPD repair"),9],
                                  y = ALL_XRcpd_enh[which(ALL_XRcpd_enh$Group=="Bottom 10% CPD - CPD repair"),9], 
                                  alternative = "two.sided", paired = F, conf.int = T)$p.value

# Correct for multiple hypothesis testing 
Enh_XRcpd_stats$Wilcoxon_RST_pval<-as.numeric(Enh_XRcpd_stats$Wilcoxon_RST_pval)

Enh_XRcpd_stats$p.adjust.BH<-p.adjust(p = Enh_XRcpd_stats$Wilcoxon_RST_pval, method = "BH")
write.csv(Enh_XRcpd_stats,file = paste0(out, "Wilcoxon RST - Enh mutR in top 10 CPD.csv"))

# Add factor levels to labels
ALL_XRcpd_enh$Group<-factor(ALL_XRcpd_enh$Group, levels = c("Top 10% CPD - no CPD repair","Top 10% CPD - CPD repair",
                                                            "Bottom 10% CPD - no CPD repair","Bottom 10% CPD - CPD repair"),
                            labels = c("Top 10% CPD -\nno CPD repair","Top 10% CPD -\nCPD repair",
                                       "Bottom 10% CPD -\nno CPD repair","Bottom 10% CPD -\nCPD repair"))

# Plot
Enh_XRcpd_All_plot<-ggplot(data = ALL_XRcpd_enh, aes(x = Group,y = Log10_mutR, fill=Group))+theme_bw()+
  geom_boxplot(show.legend = F)+ylab(yLab)+xlab(NULL)+coord_cartesian(ylim = c(0,-3.25))+
  theme(text = element_text(size = 9),axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(filename = "T10_CPD_Enh_MutRate_w_&_wo_XRcpd_repair.jpeg",
       path = out,plot = Enh_XRcpd_All_plot,device = "jpeg",width = 3,height = 4,units = "in")



# 6-4PP repair analysis - top/bottom 6-4PP-CPD regions --------------------

# Table for stats testing
RankDiff_enh_stats<-as.data.frame( matrix(nrow = 4,ncol = 3))
colnames(RankDiff_enh_stats)<-c("Group 1","Group 2","Wilcoxon_RST_pval")

# Label groups for stats test
RankDiff_enh_stats[c(1,3),1]<-"Top_10%_RankDiff-no_64PP_repair"

RankDiff_enh_stats[1,2]<-"Top_10%_6-4PP-CPD_64PP_repair"
RankDiff_enh_stats[4,1]<-"Top_10%_6-4PP-CPD_64PP_repair"

RankDiff_enh_stats[2,1]<-"Bottom_10%_6-4PP-CPD_no_64PP_repair"
RankDiff_enh_stats[3,2]<-"Bottom_10%_6-4PP-CPD_no_64PP_repair"

RankDiff_enh_stats[2,2]<-"Bottom_10%_6-4PP-CPD_64PP_repair"
RankDiff_enh_stats[4,2]<-"Bottom_10%_6-4PP-CPD_64PP_repair"


### Overlap of Enhancers in top and bottom 6-4PP-CPD regions

## Top 10%
T10_Diff_Enh_XR64<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% RankDiff")
T10_Diff_Enh_XR64$overlap<-NULL
T10_Diff_Enh_XR64$percent_overlap<-NULL
T10_Diff_Enh_XR64$Region<-NULL
T10_Diff_Enh_XR64<-as.data.frame(T10_Diff_Enh_XR64)

# Enhancers in top 10% 6-4PP-CPD regions - no 64PP repair
T10_Diff_no_XR64<-T10_Diff_Enh_XR64[which(T10_Diff_Enh_XR64$Total_XR64_strand_mean_repair==0),]
T10_Diff_no_XR64$Group<-"Top 10% 6-4PP-CPD - no 64PP repair"

# Enhancers in top 10% 6-4PP-CPD regions - with 64PP repair
T10_Diff_w_XR64<-T10_Diff_Enh_XR64[which(T10_Diff_Enh_XR64$Total_XR64_strand_mean_repair!=0),]
T10_Diff_w_XR64$Group<-"Top 10% 6-4PP-CPD - 64PP repair"


## Bottom 10%
B10_Diff_Enh_XR64<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% RankDiff")
B10_Diff_Enh_XR64$overlap<-NULL
B10_Diff_Enh_XR64$percent_overlap<-NULL
B10_Diff_Enh_XR64$Region<-NULL
B10_Diff_Enh_XR64<-as.data.frame(B10_Diff_Enh_XR64)

# Enhancers in bottom 10% 6-4PP-CPD regions - no 64PP repair
B10_Diff_no_XR64<-B10_Diff_Enh_XR64[which(B10_Diff_Enh_XR64$Total_XR64_strand_mean_repair==0),]
B10_Diff_no_XR64$Group<-"Bottom 10% 6-4PP-CPD - no 64PP repair"

# Enhancers in bottom 10% 6-4PP-CPD regions - with 64PP repair
B10_Diff_w_XR64<-B10_Diff_Enh_XR64[which(B10_Diff_Enh_XR64$Total_XR64_strand_mean_repair!=0),]
B10_Diff_w_XR64$Group<-"Bottom 10% 6-4PP-CPD - 64PP repair"


## Master table with ALL conditions
RankDiff_enh_XR64_df<-rbind(T10_Diff_no_XR64,T10_Diff_w_XR64,
                            B10_Diff_no_XR64,B10_Diff_w_XR64)
unique(RankDiff_enh_XR64_df[,10])


# Conduct Wilcoxon rank sum test
RankDiff_enh_stats[1,3]<-wilcox.test(x = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Top 10% 6-4PP-CPD - no 64PP repair"),9],
                                     y = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Top 10% 6-4PP-CPD - 64PP repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[2,3]<-wilcox.test(x = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Bottom 10% 6-4PP-CPD - no 64PP repair"),9],
                                     y = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Bottom 10% 6-4PP-CPD - 64PP repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[3,3]<-wilcox.test(x = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Top 10% 6-4PP-CPD - no 64PP repair"),9],
                                     y = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Bottom 10% 6-4PP-CPD - no 64PP repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[4,3]<-wilcox.test(x = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Top 10% 6-4PP-CPD - 64PP repair"),9],
                                     y = RankDiff_enh_XR64_df[which(RankDiff_enh_XR64_df$Group=="Bottom 10% 6-4PP-CPD - 64PP repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value

# Correct for multiple hypothesis testing 
RankDiff_enh_stats$Wilcoxon_RST_pval<-as.numeric(RankDiff_enh_stats$Wilcoxon_RST_pval)

RankDiff_enh_stats$p.adjust.BH<-p.adjust(p = RankDiff_enh_stats$Wilcoxon_RST_pval, method = "BH")
write.csv(RankDiff_enh_stats,file = paste0(out, "XR64 Wilcoxon RST - Enh mutR in top 10 6-4PP-CPD.csv"))

# Add factor levels to labels
RankDiff_enh_XR64_df$Group<-factor(RankDiff_enh_XR64_df$Group, levels = c("Top 10% 6-4PP-CPD - no 64PP repair","Top 10% 6-4PP-CPD - 64PP repair",
                                                                          "Bottom 10% 6-4PP-CPD - no 64PP repair","Bottom 10% 6-4PP-CPD - 64PP repair"),
                                   labels = c("Top 10% 6-4PP-CPD -\nno repair","Top 10% 6-4PP-CPD -\n64PP repair",
                                              "Bottom 10% 6-4PP-CPD -\nno repair","Bottom 10% 6-4PP-CPD -\n64PP repair"))

# Plot
RankDiff_All_XR64_plot<-ggplot(data = RankDiff_enh_XR64_df, aes(x = Group,y = Log10_mutR, fill=Group))+theme_bw()+
  geom_boxplot(show.legend = F)+ylab(yLab)+xlab(NULL)+coord_cartesian(ylim = c(0,-3.25))+
  theme(text = element_text(size = 9),axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(filename = "T10_RankDiff_Enh_MutRate_w_&_wo_XR64_repair.jpeg",
       path = out,plot = RankDiff_All_XR64_plot,device = "jpeg",width = 3,height = 4,units = "in")


# CPD repair analysis - top/bottom 6-4PP-CPD regions ------------------------

# Table for stats testing
RankDiff_enh_stats<-as.data.frame( matrix(nrow = 4,ncol = 3))
colnames(RankDiff_enh_stats)<-c("Group 1","Group 2","Wilcoxon_RST_pval")

# Label groups for stats test
RankDiff_enh_stats[c(1,3),1]<-"Top_10%_RankDiff-no_CPD_repair"

RankDiff_enh_stats[1,2]<-"Top_10%_RankDiff-CPD_repair"
RankDiff_enh_stats[4,1]<-"Top_10%_RankDiff-CPD_repair"

RankDiff_enh_stats[2,1]<-"Bottom_10%_RankDiff-no_CPD_repair"
RankDiff_enh_stats[3,2]<-"Bottom_10%_RankDiff-no_CPD_repair"

RankDiff_enh_stats[2,2]<-"Bottom_10%_RankDiff-CPD_repair"
RankDiff_enh_stats[4,2]<-"Bottom_10%_RankDiff-CPD_repair"
RankDiff_enh_stats


### Overlap of Enhancers in top and bottom 6-4PP-CPD regions

## Top 10%
T10_Diff_Enh_XRcpd<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% 6-4PP-CPD")
T10_Diff_Enh_XRcpd$overlap<-NULL
T10_Diff_Enh_XRcpd$percent_overlap<-NULL
T10_Diff_Enh_XRcpd$Region<-NULL
T10_Diff_Enh_XRcpd<-as.data.frame(T10_Diff_Enh_XRcpd)

# Enhancers in top 10% 6-4PP-CPD regions - no CPD repair
T10_Diff_no_XRcpd<-T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),]
T10_Diff_no_XRcpd$Group<-"Top 10% 6-4PP-CPD - no CPD repair"

# Enhancers in top 10% 6-4PP-CPD regions - with CPD repair
T10_Diff_w_XRcpd<-T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),]
T10_Diff_w_XRcpd$Group<-"Top 10% 6-4PP-CPD - CPD repair"

## Bottom 10%
B10_Diff_Enh_XRcpd<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% 6-4PP-CPD")
B10_Diff_Enh_XRcpd$overlap<-NULL
B10_Diff_Enh_XRcpd$percent_overlap<-NULL
B10_Diff_Enh_XRcpd$Region<-NULL
B10_Diff_Enh_XRcpd<-as.data.frame(B10_Diff_Enh_XRcpd)

# Enhancers in bottom 10% 6-4PP-CPD regions - no CPD repair
B10_Diff_no_XRcpd<-B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),]
B10_Diff_no_XRcpd$Group<-"Bottom 10% 6-4PP-CPD - no CPD repair"

# Enhancers in bottom 10% 6-4PP-CPD regions - with CPD repair
B10_Diff_w_XRcpd<-B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),]
B10_Diff_w_XRcpd$Group<-"Bottom 10% 6-4PP-CPD - CPD repair"

## Master table with ALL conditions
RankDiff_enh_XRcpd_df<-rbind(T10_Diff_no_XRcpd,T10_Diff_w_XRcpd,
                             B10_Diff_no_XRcpd,B10_Diff_w_XRcpd)

# Conduct Wilcoxon rank sum test
RankDiff_enh_stats[1,3]<-wilcox.test(x = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Top 10% 6-4PP-CPD - no CPD repair"),9],
                                     y = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Top 10% 6-4PP-CPD - CPD repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[2,3]<-wilcox.test(x = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Bottom 10% 6-4PP-CPD - no CPD repair"),9],
                                     y = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Bottom 10% 6-4PP-CPD - CPD repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[3,3]<-wilcox.test(x = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Top 10% 6-4PP-CPD - no CPD repair"),9],
                                     y = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Bottom 10% 6-4PP-CPD - no CPD repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value
RankDiff_enh_stats[4,3]<-wilcox.test(x = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Top 10% 6-4PP-CPD - CPD repair"),9],
                                     y = RankDiff_enh_XRcpd_df[which(RankDiff_enh_XRcpd_df$Group=="Bottom 10% 6-4PP-CPD - CPD repair"),9], 
                                     alternative = "two.sided", paired = F, conf.int = T)$p.value

# Correct for multiple hypothesis testing 
RankDiff_enh_stats$Wilcoxon_RST_pval<-as.numeric(RankDiff_enh_stats$Wilcoxon_RST_pval)

RankDiff_enh_stats$p.adjust.BH<-p.adjust(p = RankDiff_enh_stats$Wilcoxon_RST_pval, method = "BH")
write.csv(RankDiff_enh_stats,file = paste0(out, "XRcpd Wilcoxon RST - Enh mutR in top 10 6-4PP-CPD.csv"))

# Add factor levels to labels
RankDiff_enh_XRcpd_df$Group<-factor(RankDiff_enh_XRcpd_df$Group, levels = c("Top 10% 6-4PP-CPD - no CPD repair","Top 10% 6-4PP-CPD - CPD repair",
                                                                            "Bottom 10% 6-4PP-CPD - no CPD repair","Bottom 10% 6-4PP-CPD - CPD repair"),
                                    labels = c("Top 10% 6-4PP-CPD -\nno repair","Top 10% 6-4PP-CPD -\nCPD repair",
                                               "Bottom 10% 6-4PP-CPD -\nno repair","Bottom 10% 6-4PP-CPD -\nCPD repair"))

# Plot
RankDiff_All_XRcpd_plot<-ggplot(data = RankDiff_enh_XRcpd_df, aes(x = Group,y = Log10_mutR, fill=Group))+theme_bw()+
  geom_boxplot(show.legend = F)+ylab(yLab)+xlab(NULL)+coord_cartesian(ylim = c(0,-3.25))+
  theme(text = element_text(size = 9),axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(filename = "T10_RankDiff_Enh_MutRate_w_&_wo_XRcpd_repair.jpeg",
       path = out,plot = RankDiff_All_XRcpd_plot,device = "jpeg",width = 3,height = 4,units = "in")



### END ###
library(GenomicRanges)
library(ggplot2)
library(reshape2)


# Description: This analysis will compare the C>T mutation density within the top and bottom 10% regions of interest with either the whole genome
# 100kb genomic bin size


out<-"output_directory/"


# Stats function
Stats_fun<-function(Data, Ctr_group, Sample_names,ID_col,data_col) {
  Stats_Table<-as.data.frame( matrix(nrow = length(Sample_names),ncol = 5))
  Stats_Table[,1]<-Sample_names
  colnames(Stats_Table)<-c('Sample_name', 'Total_Obs_in_population','Observations_per_sample', 'Wilcoxon.RST.pval','p.adjust.BH')
  
  Stats_Table[,2]<-nrow(Data[which(Data[,ID_col]==Ctr_group),])
  
  for (i in 1:length(Sample_names)) {
    Stats_Table[i,3]<-nrow(Data[which(Data[,ID_col]==Sample_names[i]),])
    S_test<-wilcox.test(x = Data[which(Data[,ID_col]==Sample_names[i]),data_col],y = Data[which(Data[,ID_col]==Ctr_group),data_col], alternative = "two.sided", paired = F, conf.int = T)
    Stats_Table[i,4]<-S_test$p.value
  }
  Stats_Table[,5]<-p.adjust(p = Stats_Table[,4], method = "BH")
  
  return(Stats_Table)
}


# Boxplot function
Gplot_box<- function(df,Col_Name,Y_axis_var_name,Y_axis_var_col_num,
                     Group,Group_label,xLab,yLab,Margin,LEGEND) {
  All_Min<-list()
  All_Max<-list()
  
  for (i in 1:length(Group)){
    Limits<-boxplot.stats(df[which(df[colnames(df)==Col_Name]==Group[i]),Y_axis_var_col_num])
    
    Min<-Limits$stats[1]
    Max<-Limits$stats[5]
    
    Min_name<-paste0("Min_",i)
    Max_name<-paste0("Max_",i)
    
    All_Min[[Min_name]]<-Min
    All_Max[[Max_name]]<-Max
  }
  
  All_Min<-unlist(All_Min)
  All_Max<-unlist(All_Max)
  
  
  Plot<-ggplot(data = df, aes_string(x = Col_Name,y=Y_axis_var_name,fill=Col_Name))+
    geom_boxplot(show.legend = LEGEND ,outlier.shape = NA,size=0.25)+
    coord_cartesian(ylim = c(min(All_Min),max(All_Max)*Margin))+
    ylab(yLab)+xlab(xLab)+theme_bw()+
    scale_x_discrete(labels= Group_label)+
    theme(axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
          legend.title = element_blank())
  return(Plot)
  
}


## Import SBS mutation dataset for: 100kb genome. Then calculate the total C>T density by combining the normalized C>T and G>A counts into one

# Load 100kb mutation density dataset
load(file = "file_path/All_SSM_rates_100kb_hg19.RData")

# Combine both C>T and G>A to get "total C>T"
Total_CtoT<-SBS_100kb
Total_CtoT$Total_C_to_T_rate<-Total_CtoT$Normalized_C_to_T_count + Total_CtoT$Normalized_G_to_A_count
mcols(Total_CtoT)<-mcols(Total_CtoT)[25]
Total_CtoT$Region<-"Whole genome"
Total_CtoT<-subset(Total_CtoT,seqnames!="chrY")
Total_CtoT<-subset(Total_CtoT,seqnames!="chrM")


## Step 2: Load rank normalized UV lesion dataset (100kb), calculate top/bottom 10% regions
# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff, keep.extra.columns = T)

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


######
### 6-4PP

## Top 64PP
T64_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = IP64_T10)
T64_mutR$Region<-"Top 10% 64PP"
T64_mutR<-as.data.frame(T64_mutR)

## Bottom 64PP
B64_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = IP64_B10)
B64_mutR$Region<-"Bottom 10% 64PP"
B64_mutR<-as.data.frame(B64_mutR)


# Combine into one df
IP64_CtoT<-as.data.frame(Total_CtoT)
IP64_CtoT<-rbind(IP64_CtoT,T64_mutR,B64_mutR)
IP64_CtoT<-na.omit(IP64_CtoT)
IP64_CtoT$Region<-factor(IP64_CtoT$Region, levels = c("Bottom 10% 64PP","Whole genome","Top 10% 64PP"))

# Stats test
Region_names<-c("Top 10% 64PP","Bottom 10% 64PP")
IP64_CtoT_stats<-Stats_fun(Data = IP64_CtoT,Ctr_group ="Whole genome",Sample_names =Region_names,ID_col = 7,data_col = 6)

# Boxplot
IP64_CtoT_plot<-Gplot_box(df = IP64_CtoT, Col_Name = "Region",Y_axis_var_name = "Total_C_to_T_rate",Y_axis_var_col_num = 6,
                          Group = c("Bottom 10% 64PP","Whole genome","Top 10% 64PP"),Group_label = c("Bottom 10% 64PP","Whole genome","Top 10% 64PP"),
                          xLab = "",yLab = "Melanoma mutation rate\n(C>T counts per 100kb bin)",Margin = 1.5,LEGEND = F)
IP64_CtoT_plot
ggsave(filename = "IP64_MutRate_in_Top_10perc.100kb.pdf",
       path = out,plot = IP64_CtoT_plot,device = "pdf",width = 2,height = 3,units = "in")



######
## CPD

# Top CPD
Tcpd_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = CPD_T10)
Tcpd_mutR$Region<-"Top 10% CPD"
Tcpd_mutR<-as.data.frame(Tcpd_mutR)

# Bottom CPD
Bcpd_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = CPD_B10)
Bcpd_mutR$Region<-"Bottom 10% CPD"
Bcpd_mutR<-as.data.frame(Bcpd_mutR)

# Combine into one df
CPD_CtoT<-as.data.frame(Total_CtoT)
CPD_CtoT<-rbind(CPD_CtoT,Tcpd_mutR,Bcpd_mutR)
CPD_CtoT<-na.omit(CPD_CtoT)
CPD_CtoT$Region<-factor(CPD_CtoT$Region, levels = c("Bottom 10% CPD","Whole genome","Top 10% CPD"))

# Stats test
Region_names<-c("Top 10% CPD","Bottom 10% CPD")
CPD_CtoT_stats<-Stats_fun(Data = CPD_CtoT,Ctr_group ="Whole genome",Sample_names =Region_names,ID_col = 7,data_col = 6)

# Boxplot
CPD_CtoT_plot<-Gplot_box(df = CPD_CtoT, Col_Name = "Region",Y_axis_var_name = "Total_C_to_T_rate",Y_axis_var_col_num = 6,
                         Group = c("Bottom 10% CPD","Whole genome","Top 10% CPD"),Group_label = c("Bottom 10% CPD","Whole genome","Top 10% CPD"),
                         xLab = "",yLab = "Melanoma mutation rate\n(C>T counts per 100kb bin)",Margin = 1.5,LEGEND = F)
CPD_CtoT_plot
ggsave(filename = "CPD_MutRate_in_Top_10perc.100kb.pdf",
       path = out,plot = CPD_CtoT_plot,device = "pdf",width = 2,height = 3,units = "in")


######
## RankDiff

# Top RankDiff
Tdiff_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = Diff_T10)
Tdiff_mutR$Region<-"Top 10% RankDiff"
Tdiff_mutR<-as.data.frame(Tdiff_mutR)

# Bottom 64PP
Bdiff_mutR<-subsetByOverlaps(x = Total_CtoT,ranges = Diff_B10)
Bdiff_mutR$Region<-"Bottom 10% RankDiff"
Bdiff_mutR<-as.data.frame(Bdiff_mutR)

# Combine into one df
RankDiff_CtoT<-as.data.frame(Total_CtoT)
RankDiff_CtoT<-rbind(RankDiff_CtoT,Tdiff_mutR,Bdiff_mutR)
RankDiff_CtoT<-na.omit(RankDiff_CtoT)
RankDiff_CtoT$Region<-factor(RankDiff_CtoT$Region, levels = c("Bottom 10% RankDiff","Whole genome","Top 10% RankDiff"))

# Stats test
Region_names<-c("Top 10% RankDiff","Bottom 10% RankDiff")
RankDiff_CtoT_stats<-Stats_fun(Data = RankDiff_CtoT,Ctr_group ="Whole genome",Sample_names =Region_names,ID_col = 7,data_col = 6)

# Boxplot
RankDiff_CtoT_plot<-Gplot_box(df = RankDiff_CtoT, Col_Name = "Region",Y_axis_var_name = "Total_C_to_T_rate",Y_axis_var_col_num = 6,
                              Group = c("Bottom 10% RankDiff","Whole genome","Top 10% RankDiff"),Group_label = c("Bottom 10% RankDiff","Whole genome","Top 10% RankDiff"),
                              xLab = "",yLab = "Melanoma mutation rate\n(C>T counts per 100kb bin)",Margin = 1.5,LEGEND = F)
RankDiff_CtoT_plot
ggsave(filename = "RankDiff_MutRate_in_Top_10perc.100kb.pdf",
       path = out,plot = RankDiff_CtoT_plot,device = "pdf",width = 2,height = 3,units = "in")


# Save stats results:
IP64_CtoT_stats
CPD_CtoT_stats
RankDiff_CtoT_stats

All_stats<-rbind(IP64_CtoT_stats,CPD_CtoT_stats,RankDiff_CtoT_stats)

write.csv(All_stats, file = paste0(out, "Wilcox_RST_results_of_top_10perc_100kb_regions.csv"))



## END ##
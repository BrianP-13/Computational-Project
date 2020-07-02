library(GenomicRanges)
library(ggplot2)
library(IDPmisc)


## Description: This script measures the cumulative repair levels of enhancers in top regions of interest.


out<-"output_directory/"



# stats function:
Stats_fun<-function(Data, Ctr_group, Sample_names,ID_col,data_col) {
  Stats_Table<-as.data.frame( matrix(nrow = length(Sample_names),ncol = 4))
  Stats_Table[,1]<-Sample_names
  colnames(Stats_Table)<-c('Sample_name', 'Observations_per_sample', 'Wilcoxon.RST.pval','p.adjust.BH')
  
  for (i in 1:length(Sample_names)) {
    Stats_Table[i,2]<-nrow(Data[which(Data[,ID_col]==Sample_names[i]),])
    S_test<-wilcox.test(x = Data[which(Data[,ID_col]==Sample_names[i]),data_col],y = Data[which(Data[,ID_col]==Ctr_group),data_col], alternative = "two.sided", paired = F, conf.int = T)
    Stats_Table[i,3]<-S_test$p.value
  }
  Stats_Table[,4]<-p.adjust(p = Stats_Table[,3], method = "BH")
  
  return(Stats_Table)
}

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

# Modified boxplot function
Gplot_box<- function(df,Col_Name,Y_axis_var_name,Y_axis_var_col_num,
                     Group_label,xLab,yLab,Margin,LEGEND) {
  
  Plot<-ggplot(data = df, aes_string(x = Col_Name,y=Y_axis_var_name,fill=Col_Name))+
    geom_boxplot(show.legend = LEGEND ,outlier.shape = NA,size=0.25)+
    coord_cartesian(ylim = c(0,100))+
    ylab(yLab)+xlab(xLab)+theme_bw()+
    scale_x_discrete(labels= Group_label)+
    theme(axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
          legend.title = element_blank())
  return(Plot)
  
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
load(file = "file_path/Enh_Total_XR_and_mutR_dataset.RData")
Enh_XR_mutR<-subset(Enh_XR_mutR,seqnames!="chrX")

######
## 64PP repair in top 64PP regions
Enhancer_64PP_repair_table<-as.data.frame( matrix(nrow = 3,ncol = 5))
colnames(Enhancer_64PP_repair_table)<-c("Enhancer_region","Number_of_enhancers","Enhancers_with_64PP_repair","Enhancers_without_64PP_repair","Percent_enh_without_64PP_repair")
Enhancer_64PP_repair_table[,1]<-c("All enhancers","Top 10% 64PP","Bottom 10% 64PP")
Enhancer_64PP_repair_table

# Subset enhancer regions to include only 64PP repair
Enh_XR64<-Enh_XR_mutR[,1:7]
Enh_XR64$Region<-"All enhancers" # column 8

Enhancer_64PP_repair_table[1,2]<-nrow(Enh_XR64)
Enhancer_64PP_repair_table[1,3]<-nrow(Enh_XR64[which(Enh_XR64$Total_XR64_strand_mean_repair!=0),])
Enhancer_64PP_repair_table[1,4]<-nrow(Enh_XR64[which(Enh_XR64$Total_XR64_strand_mean_repair==0),])
Enhancer_64PP_repair_table[1,5]<-(Enhancer_64PP_repair_table[1,4]/Enhancer_64PP_repair_table[1,2])*100

## Subset genes overlapping in top regions
# Top 64PP
XR64_T10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_T10,Overlap = 50,Region_name = "Top 10% 64PP")
XR64_T10$overlap<-NULL
XR64_T10$percent_overlap<-NULL
XR64_T10<-as.data.frame(XR64_T10)

Enhancer_64PP_repair_table[2,2]<-nrow(XR64_T10)
Enhancer_64PP_repair_table[2,3]<-nrow(XR64_T10[which(XR64_T10$Total_XR64_strand_mean_repair!=0),])
Enhancer_64PP_repair_table[2,4]<-nrow(XR64_T10[which(XR64_T10$Total_XR64_strand_mean_repair==0),])
Enhancer_64PP_repair_table[2,5]<-(Enhancer_64PP_repair_table[2,4]/Enhancer_64PP_repair_table[2,2])*100

# Bottom 64PP
XR64_B10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_B10,Overlap = 50,Region_name = "Bottom 10% 64PP")
XR64_B10$overlap<-NULL
XR64_B10$percent_overlap<-NULL
XR64_B10<-as.data.frame(XR64_B10)

Enhancer_64PP_repair_table[3,2]<-nrow(XR64_B10)
Enhancer_64PP_repair_table[3,3]<-nrow(XR64_B10[which(XR64_B10$Total_XR64_strand_mean_repair!=0),])
Enhancer_64PP_repair_table[3,4]<-nrow(XR64_B10[which(XR64_B10$Total_XR64_strand_mean_repair==0),])
Enhancer_64PP_repair_table[3,5]<-(Enhancer_64PP_repair_table[3,4]/Enhancer_64PP_repair_table[3,2])*100

## Combine into one df
XR64_ALL<-rbind(Enh_XR64,XR64_T10,XR64_B10)
XR64_ALL<-na.omit(XR64_ALL)
XR64_ALL$Region<-factor(XR64_ALL$Region,levels = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"))

XR64_ALL<-subset(XR64_ALL,Total_XR64_strand_mean_repair!=0)

# Stats test
Region_names<-c("Top 10% 64PP","Bottom 10% 64PP")
XR64_stats<-Stats_fun(Data = XR64_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 8,data_col = 7)

# Boxplot
XR64_repair_plot<-Gplot_box(df = XR64_ALL, Col_Name = "Region",Y_axis_var_name = "Total_XR64_strand_mean_repair",Y_axis_var_col_num = 7,
                            Group_label = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"),
                            xLab = NULL,yLab = "Total 6-4PP repair",Margin = 1.25,LEGEND = F)
XR64_repair_plot
ggsave(filename = "Enhancer_64_repair_Top_10perc_64PP.pdf",
       path = out,plot = XR64_repair_plot,device = "pdf",width = 2,height = 3,units = "in")


######
## 64PP repair in top RankDiff regions
Enh_64PP_RankDiff_repair_table<-as.data.frame( matrix(nrow = 3,ncol = 5))
colnames(Enh_64PP_RankDiff_repair_table)<-c("Enhancer_region","Number_of_enhancers","Enhancers_with_64PP_repair","Enhancers_without_64PP_repair","Percent_enh_without_64PP_repair")
Enh_64PP_RankDiff_repair_table[,1]<-c("All enhancers","Top 10% RankDiff","Bottom 10% RankDiff")

Enh_64PP_RankDiff_repair_table[1,2]<-nrow(Enh_XR64)
Enh_64PP_RankDiff_repair_table[1,3]<-nrow(Enh_XR64[which(Enh_XR64$Total_XR64_strand_mean_repair!=0),])
Enh_64PP_RankDiff_repair_table[1,4]<-nrow(Enh_XR64[which(Enh_XR64$Total_XR64_strand_mean_repair==0),])
Enh_64PP_RankDiff_repair_table[1,5]<-(Enh_64PP_RankDiff_repair_table[1,4]/Enh_64PP_RankDiff_repair_table[1,2])*100


## Subset genes overlapping in top regions
# Top RankDiff
XR64_T10_Diff<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% RankDiff")
XR64_T10_Diff$overlap<-NULL
XR64_T10_Diff$percent_overlap<-NULL
XR64_T10_Diff<-as.data.frame(XR64_T10_Diff)

Enh_64PP_RankDiff_repair_table[2,2]<-nrow(XR64_T10_Diff)
Enh_64PP_RankDiff_repair_table[2,3]<-nrow(XR64_T10_Diff[which(XR64_T10_Diff$Total_XR64_strand_mean_repair!=0),])
Enh_64PP_RankDiff_repair_table[2,4]<-nrow(XR64_T10_Diff[which(XR64_T10_Diff$Total_XR64_strand_mean_repair==0),])
Enh_64PP_RankDiff_repair_table[2,5]<-(Enh_64PP_RankDiff_repair_table[2,4]/Enh_64PP_RankDiff_repair_table[2,2])*100

# Bottom 64PP
XR64_B10_Diff<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% RankDiff")
XR64_B10_Diff$overlap<-NULL
XR64_B10_Diff$percent_overlap<-NULL
XR64_B10_Diff<-as.data.frame(XR64_B10_Diff)

Enh_64PP_RankDiff_repair_table[3,2]<-nrow(XR64_B10_Diff)
Enh_64PP_RankDiff_repair_table[3,3]<-nrow(XR64_B10_Diff[which(XR64_B10_Diff$Total_XR64_strand_mean_repair!=0),])
Enh_64PP_RankDiff_repair_table[3,4]<-nrow(XR64_B10_Diff[which(XR64_B10_Diff$Total_XR64_strand_mean_repair==0),])
Enh_64PP_RankDiff_repair_table[3,5]<-(Enh_64PP_RankDiff_repair_table[3,4]/Enh_64PP_RankDiff_repair_table[3,2])*100

## Combine into one df
XR64_Diff_ALL<-rbind(Enh_XR64,XR64_T10_Diff,XR64_B10_Diff)
XR64_Diff_ALL<-na.omit(XR64_Diff_ALL)
XR64_Diff_ALL$Region<-factor(XR64_Diff_ALL$Region,levels = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"))

XR64_Diff_ALL<-subset(XR64_Diff_ALL,Total_XR64_strand_mean_repair!=0)

# Stats test
Region_names<-c("Top 10% RankDiff","Bottom 10% RankDiff")
XR64_Diff_stats<-Stats_fun(Data = XR64_Diff_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 8,data_col = 7)

# Boxplot
XR64_Diff_repair_plot<-Gplot_box(df = XR64_Diff_ALL, Col_Name = "Region",Y_axis_var_name = "Total_XR64_strand_mean_repair",Y_axis_var_col_num = 7,
                                 #Group = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),
                                 Group_label = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),
                                 xLab = NULL,yLab = "Total 6-4PP repair",Margin = 1.25,LEGEND = F)
XR64_Diff_repair_plot
ggsave(filename = "Enhancer_64_repair_Top_10perc_RankDiff.pdf",
       path = out,plot = XR64_Diff_repair_plot,device = "pdf",width = 2,height = 3,units = "in")



######
## CPD repair in top CPD regions
Enhancer_CPD_repair_table<-as.data.frame( matrix(nrow = 3,ncol = 5))
colnames(Enhancer_CPD_repair_table)<-c("Enhancer_region","Number_of_enhancers","Enhancers_with_CPD_repair","Enhancers_without_CPD_repair","Percent_enh_without_CPD_repair")
Enhancer_CPD_repair_table[,1]<-c("All enhancers","Top 10% CPD","Bottom 10% CPD")

# Subset enhancer regions to include only CPD repair
Enh_XRCPD<-Enh_XR_mutR[,c(1:6,8)]
Enh_XRCPD$Region<-"All enhancers" # column 8

Enhancer_CPD_repair_table[1,2]<-nrow(Enh_XRCPD)
Enhancer_CPD_repair_table[1,3]<-nrow(Enh_XRCPD[which(Enh_XRCPD$Total_XRCPD_strand_mean_repair!=0),])
Enhancer_CPD_repair_table[1,4]<-nrow(Enh_XRCPD[which(Enh_XRCPD$Total_XRCPD_strand_mean_repair==0),])
Enhancer_CPD_repair_table[1,5]<-(Enhancer_CPD_repair_table[1,4]/Enhancer_CPD_repair_table[1,2])*100

## Subset genes overlapping in top regions
# Top CPD
XRCPD_T10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRCPD,keep.extra.columns = T),Subject = CPD_T10,Overlap = 50,Region_name = "Top 10% CPD")
XRCPD_T10$overlap<-NULL
XRCPD_T10$percent_overlap<-NULL
XRCPD_T10<-as.data.frame(XRCPD_T10)

Enhancer_CPD_repair_table[2,2]<-nrow(XRCPD_T10)
Enhancer_CPD_repair_table[2,3]<-nrow(XRCPD_T10[which(XRCPD_T10$Total_XRCPD_strand_mean_repair!=0),])
Enhancer_CPD_repair_table[2,4]<-nrow(XRCPD_T10[which(XRCPD_T10$Total_XRCPD_strand_mean_repair==0),])
Enhancer_CPD_repair_table[2,5]<-(Enhancer_CPD_repair_table[2,4]/Enhancer_CPD_repair_table[2,2])*100

# Bottom CPD
XRCPD_B10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRCPD,keep.extra.columns = T),Subject = CPD_B10,Overlap = 50,Region_name = "Bottom 10% CPD")
XRCPD_B10$overlap<-NULL
XRCPD_B10$percent_overlap<-NULL
XRCPD_B10<-as.data.frame(XRCPD_B10)

Enhancer_CPD_repair_table[3,2]<-nrow(XRCPD_B10)
Enhancer_CPD_repair_table[3,3]<-nrow(XRCPD_B10[which(XRCPD_B10$Total_XRCPD_strand_mean_repair!=0),])
Enhancer_CPD_repair_table[3,4]<-nrow(XRCPD_B10[which(XRCPD_B10$Total_XRCPD_strand_mean_repair==0),])
Enhancer_CPD_repair_table[3,5]<-(Enhancer_CPD_repair_table[3,4]/Enhancer_CPD_repair_table[3,2])*100

## Combine into one df
XRCPD_ALL<-rbind(Enh_XRCPD,XRCPD_T10,XRCPD_B10)
XRCPD_ALL<-na.omit(XRCPD_ALL)
XRCPD_ALL$Region<-factor(XRCPD_ALL$Region,levels = c("Bottom 10% CPD","All enhancers","Top 10% CPD"))

XRCPD_ALL<-subset(XRCPD_ALL,Total_XRCPD_strand_mean_repair!=0)

# Stats test
Region_names<-c("Top 10% CPD","Bottom 10% CPD")
XRCPD_stats<-Stats_fun(Data = XRCPD_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 8,data_col = 7)
XRCPD_stats

# Boxplot
CPD_repair_plot<-Gplot_box(df = XRCPD_ALL, Col_Name = "Region",Y_axis_var_name = "Total_XRCPD_strand_mean_repair",Y_axis_var_col_num = 7,
                           #Group = c("Bottom 10% CPD","All enhancers","Top 10% CPD"),
                           Group_label = c("Bottom 10% CPD","All enhancers","Top 10% CPD"),
                           xLab = NULL,yLab = "Total CPD repair",Margin = 1.25,LEGEND = F)
CPD_repair_plot
ggsave(filename = "Enhancer_CPD_repair_Top_10perc_CPD.pdf",
       path = out,plot = CPD_repair_plot,device = "pdf",width = 2,height = 3,units = "in")


######
## CPD repair in top RankDiff regions
Enh_CPD_RankDiff_repair_table<-as.data.frame( matrix(nrow = 3,ncol = 5))
colnames(Enh_CPD_RankDiff_repair_table)<-c("Enhancer_region","Number_of_enhancers","Enhancers_with_CPD_repair","Enhancers_without_CPD_repair","Percent_enh_without_CPD_repair")
Enh_CPD_RankDiff_repair_table[,1]<-c("All enhancers","Top 10% RankDiff","Bottom 10% RankDiff")

Enh_CPD_RankDiff_repair_table[1,2]<-nrow(Enh_XRCPD)
Enh_CPD_RankDiff_repair_table[1,3]<-nrow(Enh_XRCPD[which(Enh_XRCPD$Total_XRCPD_strand_mean_repair!=0),])
Enh_CPD_RankDiff_repair_table[1,4]<-nrow(Enh_XRCPD[which(Enh_XRCPD$Total_XRCPD_strand_mean_repair==0),])
Enh_CPD_RankDiff_repair_table[1,4]<-nrow(Enh_XRCPD[which(Enh_XRCPD$Total_XRCPD_strand_mean_repair==0),])
Enh_CPD_RankDiff_repair_table[1,5]<-(Enh_CPD_RankDiff_repair_table[1,4]/Enh_CPD_RankDiff_repair_table[1,2])*100

## Subset genes overlapping in top regions
# Top RankDiff
XRCPD_T10_Diff<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRCPD,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% RankDiff")
XRCPD_T10_Diff$overlap<-NULL
XRCPD_T10_Diff$percent_overlap<-NULL
XRCPD_T10_Diff<-as.data.frame(XRCPD_T10_Diff)

Enh_CPD_RankDiff_repair_table[2,2]<-nrow(XRCPD_T10_Diff)
Enh_CPD_RankDiff_repair_table[2,3]<-nrow(XRCPD_T10_Diff[which(XRCPD_T10_Diff$Total_XRCPD_strand_mean_repair!=0),])
Enh_CPD_RankDiff_repair_table[2,4]<-nrow(XRCPD_T10_Diff[which(XRCPD_T10_Diff$Total_XRCPD_strand_mean_repair==0),])
Enh_CPD_RankDiff_repair_table[2,5]<-(Enh_CPD_RankDiff_repair_table[2,4]/Enh_CPD_RankDiff_repair_table[2,2])*100

# Bottom RankDiff
XRCPD_B10_Diff<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRCPD,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% RankDiff")
XRCPD_B10_Diff$overlap<-NULL
XRCPD_B10_Diff$percent_overlap<-NULL
XRCPD_B10_Diff<-as.data.frame(XRCPD_B10_Diff)

Enh_CPD_RankDiff_repair_table[3,2]<-nrow(XRCPD_B10_Diff)
Enh_CPD_RankDiff_repair_table[3,3]<-nrow(XRCPD_B10_Diff[which(XRCPD_B10_Diff$Total_XRCPD_strand_mean_repair!=0),])
Enh_CPD_RankDiff_repair_table[3,4]<-nrow(XRCPD_B10_Diff[which(XRCPD_B10_Diff$Total_XRCPD_strand_mean_repair==0),])
Enh_CPD_RankDiff_repair_table[3,5]<-(Enh_CPD_RankDiff_repair_table[3,4]/Enh_CPD_RankDiff_repair_table[3,2])*100

## Combine into one df
XRCPD_Diff_ALL<-rbind(Enh_XRCPD,XRCPD_T10_Diff,XRCPD_B10_Diff)
XRCPD_Diff_ALL<-na.omit(XRCPD_Diff_ALL)
XRCPD_Diff_ALL$Region<-factor(XRCPD_Diff_ALL$Region,levels = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"))

XRCPD_Diff_ALL<-subset(XRCPD_Diff_ALL,Total_XRCPD_strand_mean_repair!=0)

# Stats test
Region_names<-c("Top 10% RankDiff","Bottom 10% RankDiff")
XRCPD_Diff_stats<-Stats_fun(Data = XRCPD_Diff_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 8,data_col = 7)

# Boxplot
CPD_Diff_repair_plot<-Gplot_box(df = XRCPD_Diff_ALL, Col_Name = "Region",Y_axis_var_name = "Total_XRCPD_strand_mean_repair",Y_axis_var_col_num = 7,
                                Group_label = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),
                                xLab = NULL,yLab = "Total CPD repair",Margin = 1.25,LEGEND = F)
CPD_Diff_repair_plot
ggsave(filename = "Enhancer_CPD_repair_Top_10perc_RankDiff.pdf",
       path = out,plot = CPD_Diff_repair_plot,device = "pdf",width = 2,height = 3,units = "in")


######

## Save stats results

## 64PP repair
# Enhs with and without 64PP repair tables
Enhancer_64repair_table<-rbind(Enhancer_64PP_repair_table,Enh_64PP_RankDiff_repair_table)
Enhancer_64repair_table<-Enhancer_64repair_table[-4,]

write.csv(Enhancer_64repair_table, file = paste0(out, "Enhancer_64PP_repair_table_per_10perc_regions.csv"))

# Enh 64PP repair stats results in top regions
ALL_XR64_stats<-rbind(XR64_stats,XR64_Diff_stats)
write.csv(ALL_XR64_stats, file = paste0(out, "Wilcox_RST_results_Enh_64_repair_in_top_10perc_100kb.csv"))


## CPD repair
# Enhs with and without CPD repair tables
Enhancer_CPDrepair_table<-rbind(Enhancer_CPD_repair_table,Enh_CPD_RankDiff_repair_table)
Enhancer_CPDrepair_table<-Enhancer_CPDrepair_table[-4,]

write.csv(Enhancer_CPDrepair_table, file = paste0(out, "Enhancer_CPD_repair_table_per_10perc_regions.csv"))

# Enh CPD repair stats results in top regions
ALL_XRCPD_stats<-rbind(XRCPD_stats,XRCPD_Diff_stats)
write.csv(ALL_XRCPD_stats, file = paste0(out, "Wilcox_RST_results_Enh_CPD_repair_in_top_10perc_100kb.csv"))


## END ##
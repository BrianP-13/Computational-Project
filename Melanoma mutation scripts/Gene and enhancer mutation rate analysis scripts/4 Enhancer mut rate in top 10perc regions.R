library(GenomicRanges)
library(ggplot2)



out<-'D:/Brian/Google Drive/Stanford University/LAB/Computational Project/Publication components/Github scripts/Figure 6 scripts/temp results/'


# Overlap function
Ovrlp_fun<-function( Query, Subject,Overlap,Region_name){
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


Enhancer_mutR_table<-as.data.frame( matrix(nrow = 7,ncol = 5))
colnames(Enhancer_mutR_table)<-c("Enhancer_region","Number_of_enhancers","Enhancers_with_CtoT_mut","Enhancers_without_mut","Percent_enh_without_mut")
Enhancer_mutR_table[,1]<-c("All enhancers","Top 10% 64PP","Bottom 10% 64PP","Top 10% CPD","Bottom 10% CPD","Top 10% RankDiff","Bottom 10% RankDiff")
Enhancer_mutR_table

# Load enhancer dataset 
load(file = 'D:/Brian/Google Drive/Stanford University/LAB/Computational Project/Genomics Data/Enhancers/EnhAtlas dataset w repair and mutR/Enh_Total_XR_and_mutR_dataset.RData')
Enh_mutR<-Enh_XR_mutR
Enh_mutR$Total_XR64_strand_mean_repair<-NULL
Enh_mutR$Total_XRCPD_strand_mean_repair<-NULL
Enh_mutR$Log10_mutR<-log10(Enh_mutR$Total_C_to_T_rate)
head(Enh_mutR)

# Convert -Inf values to NA
Enh_mutR[which(Enh_mutR$Log10_mutR== "-Inf"),8]<-NA
head(Enh_mutR)

Enhancer_mutR_table[1,2]<-nrow(Enh_mutR)
Enhancer_mutR_table[1,3]<-NROW(Enh_mutR[which(Enh_mutR$Total_C_to_T_rate!=0),]) # 67800 enhancers with C>T mutations
Enhancer_mutR_table[1,4]<-NROW(Enh_mutR[which(Enh_mutR$Total_C_to_T_rate==0),]) # 17766 enhancers without C>T mutations
Enhancer_mutR_table[1,5]<-(Enhancer_mutR_table[1,4]/Enhancer_mutR_table[1,2])*100

Enh_mutR$Region<-"All enhancers"
Enh_mutR<-makeGRangesFromDataFrame(Enh_mutR,keep.extra.columns = T)
Enh_mutR

## load rank normalized UV lesion files with ranked difference column
load(file = "D:/Brian/Google Drive/Stanford University/LAB/Computational Project/DNA Lesion Analysis/Rank normalized analyses/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff, keep.extra.columns = T)
Ranked_diff

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


######
## 6-4PP

# Top 10% 64PP
Enh_T10_64<-Ovrlp_fun(Query = Enh_mutR, Subject = IP64_T10,Overlap = 50,Region_name = "Top 10% 64PP")
Enh_T10_64$overlap<-NULL
Enh_T10_64$percent_overlap<-NULL
Enh_T10_64<-as.data.frame(Enh_T10_64)
head(Enh_T10_64)

Enhancer_mutR_table[2,2]<-nrow(Enh_T10_64)
Enhancer_mutR_table[2,3]<-nrow(Enh_T10_64[Enh_T10_64$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[2,4]<-nrow(Enh_T10_64[Enh_T10_64$Total_C_to_T_rate==0,])
Enhancer_mutR_table[2,5]<-(Enhancer_mutR_table[2,4]/Enhancer_mutR_table[2,2])*100

# Bottom 10% 64PP
Enh_B10_64<-Ovrlp_fun(Query = Enh_mutR, Subject = IP64_B10,Overlap = 50,Region_name = "Bottom 10% 64PP")
Enh_B10_64$overlap<-NULL
Enh_B10_64$percent_overlap<-NULL
Enh_B10_64<-as.data.frame(Enh_B10_64)
head(Enh_B10_64)

Enhancer_mutR_table[3,2]<-nrow(Enh_B10_64)
Enhancer_mutR_table[3,3]<-nrow(Enh_B10_64[Enh_B10_64$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[3,4]<-nrow(Enh_B10_64[Enh_B10_64$Total_C_to_T_rate==0,])
Enhancer_mutR_table[3,5]<-(Enhancer_mutR_table[3,4]/Enhancer_mutR_table[3,2])*100

# Combine top and bottom regions into one df
Enh_64_ALL<-rbind(as.data.frame(Enh_mutR),Enh_T10_64,Enh_B10_64)
Enh_64_ALL$Region<-factor(Enh_64_ALL$Region, levels = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"))
head(Enh_64_ALL)

# Stats test
Region_names<-c("Top 10% 64PP","Bottom 10% 64PP")

# Stats test using Log10 mut rate
colnames(Enh_64_ALL)
Enh_64_ALL<-na.omit(Enh_64_ALL)

Enh_64_stats<-Stats_fun(Data = Enh_64_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 9,data_col = 8)
Enh_64_stats

# Boxplot using non-transformed mut rate
Enh_64_plot<-Gplot_box(df = Enh_64_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                       Group = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"),Group_label = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"),
                       xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
Enh_64_plot
ggsave(filename = "Enhancer_mutRate_in_Top10_64PP.pdf",
       path = out,plot = Enh_64_plot,device = "pdf",width = 2,height = 3,units = "in")

# Boxplot with no labels
Enh_64_plot2<-Gplot_box(df = Enh_64_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                       Group = c("Bottom 10% 64PP","All enhancers","Top 10% 64PP"),Group_label = NULL,
                       xLab = NULL,yLab = NULL,Margin = 0.35,LEGEND = F)+scale_y_continuous(labels = NULL, minor_breaks = NULL)+
  theme_bw()+theme(line = element_line(size = 0.25), 
                   panel.border = element_rect(fill = NA,size = 0.25),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank())
Enh_64_plot2
ggsave(filename = "Enhancer_mutRate_in_Top10_64PP.NO.LABS.pdf",
       path = out,plot = Enh_64_plot2,device = "pdf",width = 1,height = 1.5,units = "in")


##

ggsave(filename = ".pdf",
       path = out,plot = ,device = "pdf",width = 2,height = 3,units = "in")

theme_bw()+theme(line = element_line(size = 0.25), 
                 panel.border = element_rect(fill = NA,size = 0.25),
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_blank())

ggsave(filename = ".NO.LABS.pdf",
       path = out,plot = ,device = "pdf",width = 1,height = 1.5,units = "in")





######
## CPD

# Top 10% CPD
Enh_T10_CPD<-Ovrlp_fun(Query = Enh_mutR, Subject = CPD_T10,Overlap = 50,Region_name = "Top 10% CPD")
Enh_T10_CPD$overlap<-NULL
Enh_T10_CPD$percent_overlap<-NULL
Enh_T10_CPD<-as.data.frame(Enh_T10_CPD)

Enhancer_mutR_table[4,2]<-nrow(Enh_T10_CPD)
Enhancer_mutR_table[4,3]<-nrow(Enh_T10_CPD[Enh_T10_CPD$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[4,4]<-nrow(Enh_T10_CPD[Enh_T10_CPD$Total_C_to_T_rate==0,])
Enhancer_mutR_table[4,5]<-(Enhancer_mutR_table[4,4]/Enhancer_mutR_table[4,2])*100

# Bottom 10% CPD
Enh_B10_CPD<-Ovrlp_fun(Query = Enh_mutR, Subject = CPD_B10,Overlap = 50,Region_name = "Bottom 10% CPD")
Enh_B10_CPD$overlap<-NULL
Enh_B10_CPD$percent_overlap<-NULL
Enh_B10_CPD<-as.data.frame(Enh_B10_CPD)
head(Enh_B10_CPD)

Enhancer_mutR_table[5,2]<-nrow(Enh_B10_CPD)
Enhancer_mutR_table[5,3]<-nrow(Enh_B10_CPD[Enh_B10_CPD$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[5,4]<-nrow(Enh_B10_CPD[Enh_B10_CPD$Total_C_to_T_rate==0,])
Enhancer_mutR_table[5,5]<-(Enhancer_mutR_table[5,4]/Enhancer_mutR_table[5,2])*100


# Combine top and bottom regions into one df
Enh_CPD_ALL<-rbind(as.data.frame(Enh_mutR),Enh_T10_CPD,Enh_B10_CPD)
Enh_CPD_ALL$Region<-factor(Enh_CPD_ALL$Region, levels = c("Bottom 10% CPD","All enhancers","Top 10% CPD"))
head(Enh_CPD_ALL)


# Stats test using Log10 mut rate
Region_names<-c("Top 10% CPD","Bottom 10% CPD")
Enh_CPD_ALL<-na.omit(Enh_CPD_ALL)

Enh_CPD_stats<-Stats_fun(Data = Enh_CPD_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 9,data_col = 8)
Enh_CPD_stats

# Boxplot using non-transformed mut rate
Enh_CPD_plot<-Gplot_box(df = Enh_CPD_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                       Group = c("Bottom 10% CPD","All enhancers","Top 10% CPD"),Group_label = c("Bottom 10% CPD","All enhancers","Top 10% CPD"),
                       xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
Enh_CPD_plot
ggsave(filename = "Enhancer_mutRate_in_Top10_CPD.pdf",
       path = out,plot = Enh_CPD_plot,device = "pdf",width = 2,height = 3,units = "in")

# Boxplot with no labels
Enh_CPD_plot2<-Gplot_box(df = Enh_CPD_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                        Group = c("Bottom 10% CPD","All enhancers","Top 10% CPD"),Group_label = NULL,
                        xLab = NULL,yLab = NULL,Margin = 0.35,LEGEND = F)+scale_y_continuous(labels = NULL,minor_breaks = NULL)+
  theme_bw()+theme(line = element_line(size = 0.25), 
                   panel.border = element_rect(fill = NA,size = 0.25),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank())
Enh_CPD_plot2
ggsave(filename = "Enhancer_mutRate_in_Top10_CPD.NO.LABS.pdf",
       path = out,plot = Enh_CPD_plot2,device = "pdf",width = 1,height = 1.5,units = "in")


######
## RankDiff

# Top 10% RankDiff
Enh_T10_Diff<-Ovrlp_fun(Query = Enh_mutR, Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% RankDiff")
Enh_T10_Diff$overlap<-NULL
Enh_T10_Diff$percent_overlap<-NULL
Enh_T10_Diff<-as.data.frame(Enh_T10_Diff)
head(Enh_T10_Diff)

Enhancer_mutR_table[6,2]<-nrow(Enh_T10_Diff)
Enhancer_mutR_table[6,3]<-nrow(Enh_T10_Diff[Enh_T10_Diff$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[6,4]<-nrow(Enh_T10_Diff[Enh_T10_Diff$Total_C_to_T_rate==0,])
Enhancer_mutR_table[6,5]<-(Enhancer_mutR_table[6,4]/Enhancer_mutR_table[6,2])*100

# Bottom 10% RankDiff
Enh_B10_Diff<-Ovrlp_fun(Query = Enh_mutR, Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% RankDiff")
Enh_B10_Diff$overlap<-NULL
Enh_B10_Diff$percent_overlap<-NULL
Enh_B10_Diff<-as.data.frame(Enh_B10_Diff)
head(Enh_B10_Diff)

Enhancer_mutR_table[7,2]<-nrow(Enh_B10_Diff)
Enhancer_mutR_table[7,3]<-nrow(Enh_B10_Diff[Enh_B10_Diff$Total_C_to_T_rate!=0,])
Enhancer_mutR_table[7,4]<-nrow(Enh_B10_Diff[Enh_B10_Diff$Total_C_to_T_rate==0,])
Enhancer_mutR_table[7,5]<-(Enhancer_mutR_table[7,4]/Enhancer_mutR_table[7,2])*100

# Combine top and bottom regions into one df
Enh_Diff_ALL<-rbind(as.data.frame(Enh_mutR),Enh_T10_Diff,Enh_B10_Diff)
Enh_Diff_ALL$Region<-factor(Enh_Diff_ALL$Region, levels = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"))
head(Enh_Diff_ALL)

# Stats test
Region_names<-c("Top 10% RankDiff","Bottom 10% RankDiff")
Enh_Diff_ALL<-na.omit(Enh_Diff_ALL)

# Stats test using Log10 mut rate
Enh_Diff_stats<-Stats_fun(Data = Enh_Diff_ALL,Ctr_group ="All enhancers",Sample_names =Region_names,ID_col = 9,data_col = 8)
Enh_Diff_stats

# Boxplot using non-transformed mut rate
Enh_Diff_plot<-Gplot_box(df = Enh_Diff_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                       Group = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),Group_label = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),
                       xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
Enh_Diff_plot
ggsave(filename = "Enhancer_mutRate_in_Top10_RankDiff.pdf",
       path = out,plot = Enh_Diff_plot,device = "pdf",width = 2,height = 3,units = "in")

# Boxplot with no labels
Enh_Diff_plot2<-Gplot_box(df = Enh_Diff_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mutR",Y_axis_var_col_num = 8,
                         Group = c("Bottom 10% RankDiff","All enhancers","Top 10% RankDiff"),Group_label = NULL,
                         xLab = NULL,yLab = NULL,Margin = 0.35,LEGEND = F)+scale_y_continuous(labels = NULL,minor_breaks = NULL)+
  theme_bw()+theme(line = element_line(size = 0.25), 
                   panel.border = element_rect(fill = NA,size = 0.25),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank())
Enh_Diff_plot2
ggsave(filename = "Enhancer_mutRate_in_Top10_RankDiff.NO.LABS.pdf",
       path = out,plot = Enh_Diff_plot2,device = "pdf",width = 1,height = 1.5,units = "in")


######

Enhancer_mutR_table

write.csv(Enhancer_mutR_table, file = paste0(out,"Num_of_mutated_enhancers_in_Top10_100kb_regions.csv"))


Enh_All_stats<-rbind(Enh_64_stats,Enh_CPD_stats,Enh_Diff_stats)
Enh_All_stats

write.csv(Enh_All_stats, file = paste0(out,"Wilcoxon_RST_of_Enh_in_Top10_100kb_regions.csv"))




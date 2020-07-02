library(GenomicRanges)
library(ggplot2)
library(IDPmisc)


## Description: 


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
  hits<-findOverlaps(query = Query,subject =  Subject,type = "within",ignore.strand=T,select = "all")
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


Gene_mutR_table<-as.data.frame( matrix(nrow = 7,ncol = 4))
colnames(Gene_mutR_table)<-c("Gene_region","Number_of_genes","Genes_with_CtoT_mut","Genes_without_mut")
Gene_mutR_table[,1]<-c("All genes","Top 10% 64PP","Bottom 10% 64PP","Top 10% CPD","Bottom 10% CPD","Top 10% RankDiff","Bottom 10% RankDiff")



# Import rank normalized UV dmg dataset and calculate top and bottom 2SD regions
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff,keep.extra.columns = T)

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


######
## 64PP

# Load 64PP dataset with dmg and repair
load(file = "file_path/Gene_list_64PP_dmg_repair_mutR.RData")
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,seqnames!="chrM")
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,seqnames!="chrY")

# Subset to only include UV damage and mutation rate
IP64_mutR<-IP64_dmg_repair_mutR
IP64_mutR$Total_XR64_strand_mean_repair<-NULL
IP64_mutR$Region<-"All genes"
IP64_mutR<-as.data.frame(IP64_mutR)

Gene_mutR_table[1,2]<-nrow(IP64_mutR)
Gene_mutR_table[1,3]<-nrow(IP64_mutR[which(IP64_mutR$Total_C_to_T_rate!=0),])
Gene_mutR_table[1,4]<-nrow(IP64_mutR[which(IP64_mutR$Total_C_to_T_rate==0),])

# Subset genes overlapping in top regions
IP64_T10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(IP64_mutR,keep.extra.columns = T),Subject = IP64_T10,Overlap = 50,Region_name = "Top 10% 64PP")
IP64_T10$overlap<-NULL
IP64_T10$percent_overlap<-NULL
IP64_T10<-as.data.frame(IP64_T10)

Gene_mutR_table[2,2]<-nrow(IP64_T10)
Gene_mutR_table[2,3]<-nrow(IP64_T10[which(IP64_T10$Total_C_to_T_rate!=0),])
Gene_mutR_table[2,4]<-nrow(IP64_T10[which(IP64_T10$Total_C_to_T_rate==0),])

IP64_B10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(IP64_mutR,keep.extra.columns = T),Subject = IP64_B10,Overlap = 50,Region_name = "Bottom 10% 64PP")
IP64_B10$overlap<-NULL
IP64_B10$percent_overlap<-NULL
IP64_B10<-as.data.frame(IP64_B10)

Gene_mutR_table[3,2]<-nrow(IP64_B10)
Gene_mutR_table[3,3]<-nrow(IP64_B10[which(IP64_B10$Total_C_to_T_rate!=0),])
Gene_mutR_table[3,4]<-nrow(IP64_B10[which(IP64_B10$Total_C_to_T_rate==0),])

# Combine into one df
IP64_mutR_ALL<-rbind(IP64_mutR,IP64_T10,IP64_B10)
IP64_mutR_ALL<-na.omit(IP64_mutR_ALL)
IP64_mutR_ALL<-subset(IP64_mutR_ALL,Log10_mut_rate!="-Inf")
IP64_mutR_ALL$Region<-factor(IP64_mutR_ALL$Region,levels = c("Bottom 10% 64PP","All genes","Top 10% 64PP"))

# Stats test
Region_names<-c("Top 10% 64PP","Bottom 10% 64PP")
IP64_mutR_stats<-Stats_fun(Data = IP64_mutR_ALL,Ctr_group ="All genes",Sample_names =Region_names,ID_col = 12,data_col = 11)

# Boxplot
IP64_mutR_plot<-Gplot_box(df = IP64_mutR_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mut_rate",Y_axis_var_col_num = 11,
                          Group = c("Bottom 10% 64PP","All genes","Top 10% 64PP"),Group_label = c("Bottom 10% 64PP","All genes","Top 10% 64PP"),
                          xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
IP64_mutR_plot
ggsave(filename = "Gene_MutRate_Top_10perc_64PP_100kb.pdf",
       path = out,plot = IP64_mutR_plot,device = "pdf",width = 2,height = 3,units = "in")



######
## CPD

# Load CPD dataset with dmg and repair
load(file = "file_path/Gene_list_CPD_dmg_repair_mutR.RData")
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,seqnames!="chrM")
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,seqnames!="chrY")

# Subset to only include UV damage and mutation rate
CPD_mutR<-CPD_dmg_repair_mutR
CPD_mutR$Total_XRCPD_strand_mean_repair<-NULL
CPD_mutR$Region<-"All genes"
CPD_mutR<-as.data.frame(CPD_mutR)

# Subset genes overlapping in top regions
CPD_T10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(CPD_mutR,keep.extra.columns = T),Subject = CPD_T10,Overlap = 50,Region_name = "Top 10% CPD")
CPD_T10$overlap<-NULL
CPD_T10$percent_overlap<-NULL
CPD_T10<-as.data.frame(CPD_T10)

Gene_mutR_table[4,2]<-nrow(CPD_T10)
Gene_mutR_table[4,3]<-nrow(CPD_T10[which(CPD_T10$Total_C_to_T_rate!=0),])
Gene_mutR_table[4,4]<-nrow(CPD_T10[which(CPD_T10$Total_C_to_T_rate==0),])

CPD_B10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(CPD_mutR,keep.extra.columns = T),Subject = CPD_B10,Overlap = 50,Region_name = "Bottom 10% CPD")
CPD_B10$overlap<-NULL
CPD_B10$percent_overlap<-NULL
CPD_B10<-as.data.frame(CPD_B10)

Gene_mutR_table[5,2]<-nrow(CPD_B10)
Gene_mutR_table[5,3]<-nrow(CPD_B10[which(CPD_B10$Total_C_to_T_rate!=0),])
Gene_mutR_table[5,4]<-nrow(CPD_B10[which(CPD_B10$Total_C_to_T_rate==0),])

# Combine into one df
CPD_mutR_ALL<-rbind(CPD_mutR,CPD_T10,CPD_B10)
CPD_mutR_ALL<-na.omit(CPD_mutR_ALL)
CPD_mutR_ALL<-subset(CPD_mutR_ALL,Log10_mut_rate!="-Inf")
CPD_mutR_ALL$Region<-factor(CPD_mutR_ALL$Region,levels = c("Bottom 10% CPD","All genes","Top 10% CPD"))

# Stats test
Region_names<-c("Top 10% CPD","Bottom 10% CPD")
CPD_mutR_stats<-Stats_fun(Data = CPD_mutR_ALL,Ctr_group ="All genes",Sample_names =Region_names,ID_col = 12,data_col = 11)

# Boxplot
CPD_mutR_plot<-Gplot_box(df = CPD_mutR_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mut_rate",Y_axis_var_col_num = 11,
                          Group = c("Bottom 10% CPD","All genes","Top 10% CPD"),Group_label = c("Bottom 10% CPD","All genes","Top 10% CPD"),
                          xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
CPD_mutR_plot
ggsave(filename = "Gene_MutRate_Top_10perc_CPD_100kb.pdf",
       path = out,plot = CPD_mutR_plot,device = "pdf",width = 2,height = 3,units = "in")


######
## RankDiff

# Subset genes overlapping in top regions
Diff_T10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(IP64_mutR,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% RankDiff")
Diff_T10$overlap<-NULL
Diff_T10$percent_overlap<-NULL
Diff_T10<-as.data.frame(Diff_T10)

Gene_mutR_table[6,2]<-nrow(Diff_T10)
Gene_mutR_table[6,3]<-nrow(Diff_T10[which(Diff_T10$Total_C_to_T_rate!=0),])
Gene_mutR_table[6,4]<-nrow(Diff_T10[which(Diff_T10$Total_C_to_T_rate==0),])

Diff_B10<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(IP64_mutR,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% RankDiff")
Diff_B10$overlap<-NULL
Diff_B10$percent_overlap<-NULL
Diff_B10<-as.data.frame(Diff_B10)

Gene_mutR_table[7,2]<-nrow(Diff_B10)
Gene_mutR_table[7,3]<-nrow(Diff_B10[which(Diff_B10$Total_C_to_T_rate!=0),])
Gene_mutR_table[7,4]<-nrow(Diff_B10[which(Diff_B10$Total_C_to_T_rate==0),])

# Combine into one df
Diff_mutR_ALL<-rbind(IP64_mutR,Diff_T10,Diff_B10)
Diff_mutR_ALL<-na.omit(Diff_mutR_ALL)
Diff_mutR_ALL<-subset(Diff_mutR_ALL,Log10_mut_rate!="-Inf")
Diff_mutR_ALL$Region<-factor(Diff_mutR_ALL$Region,levels = c("Bottom 10% RankDiff","All genes","Top 10% RankDiff"))

# Stats test
Region_names<-c("Top 10% RankDiff","Bottom 10% RankDiff")
Diff_mutR_stats<-Stats_fun(Data = Diff_mutR_ALL,Ctr_group ="All genes",Sample_names =Region_names,ID_col = 12,data_col = 11)

# Boxplot
Diff_mutR_plot<-Gplot_box(df = Diff_mutR_ALL, Col_Name = "Region",Y_axis_var_name = "Log10_mut_rate",Y_axis_var_col_num = 11,
                          Group = c("Bottom 10% RankDiff","All genes","Top 10% RankDiff"),Group_label = c("Bottom 10% RankDiff","All genes","Top 10% RankDiff"),
                          xLab = NULL,yLab = "Log10(total C>T mutation rate)",Margin = 0.35,LEGEND = F)
Diff_mutR_plot
ggsave(filename = "Gene_MutRate_Top_10perc_RankDiff_100kb.pdf",
       path = out,plot = Diff_mutR_plot,device = "pdf",width = 2,height = 3,units = "in")



# Save stats results:
IP64_mutR_stats
CPD_mutR_stats
Diff_mutR_stats

All_stats<-rbind(IP64_mutR_stats,CPD_mutR_stats,Diff_mutR_stats)

write.csv(All_stats, file = paste0(out, "Wilcox_RST_results_gene_mutR_top_10perc_100kb.csv"))


## END ##
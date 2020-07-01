library(ggplot2)

## Description: This script is for visualizing di-pyrimidine frequency distribution (boxplot)


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

Gplot_box<- function(df,Col_Name,Y_axis_var_name,Y_axis_var_col_num,
                     Group,Group_label,xLab,yLab,Margin,LEGEND) {
  # df = dataframe
  # Col_Name              = Character name of column with numerical data (experimental variable, continuous)
  # Y_axis_var_name       = Character. Name of column with sample names. 
  # Y_axis_var_col_num    = Column number with numerical values
  # Group                 = Character string of names/IDs for the samples to plot (will go on x-axis, e.g. c("control/WT","experimental/mutant1"))
  # Group_label           = Name of samples you want to appear on plot
  # xLab                  = X axis label
  # yLab                  = Y axis label
  # Margin                = Need addition space on y-axis? Value is numerical. "1" equals normal scaling.
  # LEGEND                = Want legend included in plot? (Logical, TRUE or FALSE)
  
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


######
## Measure di-Py frequencies within each replicate file ###
Path<-'path_to_DiPy_freq_files/'

# Input
Input<-read.csv(paste0(Path, "/64PP_pooled_Input_diPy_FREQ_FASTQ.csv", sep=""))
Input<-Input[,c(17,15,9,7)]
Input$Group<-"Input"

# IP 64PP
IP100<-read.csv(paste0(Path, "/64PP_pooled_IP100_diPy_FREQ_FASTQ.csv", sep=""))
IP100<-IP100[,c(17,15,9,7)]
IP100$Group<-"IP_64PP"

# Combine into one dataframe
DiNt_freq<-rbind(Input,IP100)
DiNt_freq$Group<-factor(DiNt_freq$Group,levels = c("Input","IP_64PP"))

## Stats test:
DiPy<-c("TT","TC"<"CT","CC")

# TT
TT_stats<-Stats_fun(Data = DiNt_freq,Ctr_group = "Input",Sample_names = c("Input","IP_64PP"),ID_col = 5,data_col = 1)
TT_stats[,1]<-paste0(TT_stats[,1],"_TT_freq")

# TC
TC_stats<-Stats_fun(Data = DiNt_freq,Ctr_group = "Input",Sample_names = c("Input","IP_64PP"),ID_col = 5,data_col = 2)
TC_stats[,1]<-paste0(TC_stats[,1],"_TC_freq")

# CT
CT_stats<-Stats_fun(Data = DiNt_freq,Ctr_group = "Input",Sample_names = c("Input","IP_64PP"),ID_col = 5,data_col = 3)
CT_stats[,1]<-paste0(CT_stats[,1],"_CT_freq")

# CC
CC_stats<-Stats_fun(Data = DiNt_freq,Ctr_group = "Input",Sample_names = c("Input","IP_64PP"),ID_col = 5,data_col = 4)
CC_stats[,1]<-paste0(CC_stats[,1],"_CC_freq")

Stats_table<-rbind(TT_stats,TC_stats,CT_stats,CC_stats)
Stats_table<-Stats_table[c(2,4,6,8),]

write.csv(Stats_table, file = paste0(out,"WRST_DiNt_freq_Input_vs_IP64PP100.csv"))



## Boxplots

# TT plots, with labels
TT_plot<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "TT",Y_axis_var_col_num = 1,
                   Group = c("Input","IP_64PP"),Group_label = c("Input","IP 6-4PP"),
                   xLab = NULL,yLab = "TT frequency",Margin = 1.125,LEGEND = F)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
TT_plot
ggsave(filename = "TT_freq_Input_vs_IP64PP100.pdf",
       path = out,plot = TT_plot,device = "pdf",width = 1,height = 1.5,units = "in")

# TT plots, NO labels
TT_plot2<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "TT",Y_axis_var_col_num = 1,
                   Group = c("Input","IP_64PP"),Group_label = NULL,
                   xLab = NULL,yLab = NULL,Margin = 1.125,LEGEND = F)+scale_y_continuous(labels = NULL, minor_breaks = NULL)+
  theme_bw()+ theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank())
TT_plot2
ggsave(filename = "TT_freq_Input_vs_IP64PP100.NO.LABS.pdf",
       path = out,plot = TT_plot2,device = "pdf",width = 0.75,height = 1.5,units = "in")

#
# TC plots, with labels
TC_plot<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "TC",Y_axis_var_col_num = 2,
                   Group = c("Input","IP_64PP"),Group_label = c("Input","IP 6-4PP"),
                   xLab = "",yLab = "TC frequency",Margin = 1.125,LEGEND = F)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
TC_plot
ggsave(filename = "TC_freq_Input_vs_IP64PP100.pdf",
       path = out,plot = TC_plot,device = "pdf",width = 1,height = 1.5,units = "in")

# TC plots, NO labels
TC_plot2<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "TC",Y_axis_var_col_num = 2,
                    Group = c("Input","IP_64PP"),Group_label = NULL,
                    xLab = NULL,yLab = NULL,Margin = 1.125,LEGEND = F)+scale_y_continuous(labels = NULL, minor_breaks = NULL)+
  theme_bw()+ theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank())
TC_plot2
ggsave(filename = "TC_freq_Input_vs_IP64PP100.NO.LABS.pdf",
       path = out,plot = TC_plot2,device = "pdf",width = 0.75,height = 1.5,units = "in")

#
# CT plots, with labels
CT_plot<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "CT",Y_axis_var_col_num = 3,
                   Group = c("Input","IP_64PP"),Group_label = c("Input","IP 6-4PP"),
                   xLab = "",yLab = "CT frequency",Margin = 1.05,LEGEND = F)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
CT_plot
ggsave(filename = "CT_freq_Input_vs_IP64PP100.pdf",
       path = out,plot = CT_plot,device = "pdf",width = 1,height = 1.5,units = "in")

# CT plots, NO labels
CT_plot2<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "CT",Y_axis_var_col_num = 3,
                    Group = c("Input","IP_64PP"),Group_label = NULL,
                    xLab = NULL,yLab = NULL,Margin = 1.05,LEGEND = F)+scale_y_continuous(labels = NULL, minor_breaks = NULL)+
  theme_bw()+ theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank())
CT_plot2
ggsave(filename = "CT_freq_Input_vs_IP64PP100.NO.LABS.pdf",
       path = out,plot = CT_plot2,device = "pdf",width = 0.75,height = 1.5,units = "in")

#
# CC plots, with labels
CC_plot<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "CC",Y_axis_var_col_num = 4,
                   Group = c("Input","IP_64PP"),Group_label = c("Input","IP 6-4PP"),
                   xLab = "",yLab = "CC frequency",Margin = 1.2,LEGEND = F)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
CC_plot
ggsave(filename = "CC_freq_Input_vs_IP64PP100.pdf",
       path = out,plot = CC_plot,device = "pdf",width = 1,height = 1.5,units = "in")

# CC plots, NO labels
CC_plot2<-Gplot_box(df = DiNt_freq,Col_Name = "Group",Y_axis_var_name = "CC",Y_axis_var_col_num = 4,
                    Group = c("Input","IP_64PP"),Group_label = NULL,
                    xLab = NULL,yLab = NULL,Margin = 1.2,LEGEND = F)+scale_y_continuous(labels = NULL, minor_breaks = NULL)+
  theme_bw()+ theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank())
CC_plot2
ggsave(filename = "CC_freq_Input_vs_IP64PP100.NO.LABS.pdf",
       path = out,plot = CC_plot2,device = "pdf",width = 0.75,height = 1.5,units = "in")


# END


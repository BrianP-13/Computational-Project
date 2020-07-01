library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(IDPmisc)


## Description: This script measures the di-pyrimidine frequency within each chromatin state and contains additional code for generating boxplots to visualize results across all 15 states

out<-"output_directory/"


Gplot_box<- function(df,Col_Name,Y_axis_var_name,Y_axis_var_col_num,
                     Group,Group_label,xLab,yLab,Fill_color,Margin,LEGEND) {
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
    geom_boxplot(fill= Fill_color ,show.legend = LEGEND ,outlier.shape = NA,size=0.25)+
    coord_cartesian(ylim = c(min(All_Min),max(All_Max)*Margin))+
    ylab(yLab)+xlab(xLab)+theme_bw()+
    scale_x_discrete(labels= Group_label)+
    theme(axis.title.y = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size=8,angle = 45,hjust = 1,vjust = 1),
          legend.title = element_blank())
  return(Plot)
  
}


## Get DNA sequence data for hg19 
BSgenome<-BSgenome.Hsapiens.UCSC.hg19
hg19<-seqinfo(BSgenome)


## Load IMR90 15 core chromatin state dataset
core15<-file.path("file_path/E017_15_coreMarks_dense.bed.bgz")
df_core15<-import.bed(core15)
mcols(df_core15)$score<-NULL
mcols(df_core15)$itemRgb<-NULL
mcols(df_core15)$thick<-NULL
df_core15<-df_core15[which(df_core15@seqnames != "chrY"),]
df_core15<-df_core15[which(df_core15@seqnames != "chrM"),]
df_core15<-sortSeqlevels(df_core15)
df_core15<-sort(df_core15)

# chromatin state names
core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_names2<-c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies")


## Calculate diNt freq for all chromatin states
DNAseqs<-getSeq(BSgenome, df_core15)

# Get all di Nt freq
Freq<-dinucleotideFrequency(DNAseqs, as.prob = T)

# Add diNt data to chrom state dataset
mcols(df_core15)<-cbind(mcols(df_core15), as.data.frame(Freq))
save(df_core15, file = paste0(out, "IMR90_chrom_states_with_DiPy_freq.RData")) # Good spot to save results and stop 


######
## Load Chromatin state file with Di-Py freq data
load(file = "file_path/IMR90_chrom_states_with_DiPy_freq.RData")
df_core15

core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_names2<-c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies")

DiNts<-c("TT","TC","CT","CC")
DiNt_chrom15<-df_core15
DiNt_chrom15<-as.data.frame(mcols(DiNt_chrom15)[c(17,15,9,7,1)])
colnames(DiNt_chrom15)[5]<-"Region"
DiNt_chrom15$Region<-"Whole genome"

for (i in 1:15) {
  DiNt_chrom_state<-subset(df_core15, name==core15_names[i])
  DiNt_chrom_state<-as.data.frame(mcols(DiNt_chrom_state)[c(17,15,9,7,1)])
  colnames(DiNt_chrom_state)[5]<-"Region"
  
  DiNt_chrom15<-rbind(DiNt_chrom15,DiNt_chrom_state)
}


## Plots using RankDiff order
XLAB<-"IMR90 15 core chromatin states"

## Factor levels using numerical order of chromatin states
Chrom_order_2<-c("15_Quies","5_TxWk","9_Het","14_ReprPCWk","Whole genome",
  "4_Tx","8_ZNF/Rpts","7_Enh","6_EnhG","1_TssA","3_TxFlnk",
  "13_ReprPC","2_TssAFlnk","10_TssBiv","12_EnhBiv","11_BivFlnk")

Chrom_labels_2<-c("Quies","TxWk","Het","ReprPCWk","Whole genome",
                  "Tx","ZNF/Rpts","Enh","EnhG","TssA","TxFlnk",
                  "ReprPC","TssAFlnk","TssBiv","EnhBiv","BivFlnk")
DiNt_chrom15$Region<-factor(DiNt_chrom15$Region,levels = Chrom_order_2)


## TT freq plots
TT_plot<-Gplot_box(df = DiNt_chrom15,Col_Name = "Region",
                   Y_axis_var_name = "TT",Y_axis_var_col_num = 1,
                   Group = Chrom_order_2,Group_label = Chrom_labels_2,
                   xLab = XLAB,yLab = "TT frequency",Fill_color = "Indian Red",Margin = 1,LEGEND = F)
TT_plot
ggsave(filename = "TT_freq_in_chrom_states.RankDiff_order.pdf",
       path = out,plot = TT_plot,device = "pdf",width = 3,height = 1.5,units = "in")

#
## TC freq plots
TC_plot<-Gplot_box(df = DiNt_chrom15,Col_Name = "Region",
                   Y_axis_var_name = "TC",Y_axis_var_col_num = 2,
                   Group = Chrom_order_2,Group_label = Chrom_labels_2,
                   xLab = XLAB,yLab = "TC frequency",Fill_color = "Pale Turquoise",Margin = 1,LEGEND = F)
TC_plot
ggsave(filename = "TC_freq_in_chrom_states.RankDiff_order.pdf",
       path = out,plot = TC_plot,device = "pdf",width = 3,height = 1.5,units = "in")

#
## CT freq plots
CT_plot<-Gplot_box(df = DiNt_chrom15,Col_Name = "Region",
                   Y_axis_var_name = "CT",Y_axis_var_col_num = 3,
                   Group = Chrom_order_2,Group_label = Chrom_labels_2,
                   xLab = XLAB,yLab = "CT frequency",Fill_color = "chocolate1",Margin = 1,LEGEND = F)
CT_plot
ggsave(filename = "CT_freq_in_chrom_states.RankDiff_order.pdf",
       path = out,plot = CT_plot,device = "pdf",width = 3,height = 1.5,units = "in")


#
## CC freq plots
CC_plot<-Gplot_box(df = DiNt_chrom15,Col_Name = "Region",
                   Y_axis_var_name = "CC",Y_axis_var_col_num = 4,
                   Group = Chrom_order_2,Group_label = Chrom_labels_2,
                   xLab = XLAB,yLab = "CC frequency",Fill_color = "darkolivegreen3",Margin = 1,LEGEND = F)
CC_plot
ggsave(filename = "CC_freq_in_chrom_states.RankDiff_order.pdf",
       path = out,plot = CC_plot,device = "pdf",width = 3,height = 1.5,units = "in")



## END ##
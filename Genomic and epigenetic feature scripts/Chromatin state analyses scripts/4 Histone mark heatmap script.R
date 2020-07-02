library(ggplot2)
library(reshape2)


## Description: Script for generating heatmap of histone mark fold change levels in IMR90 15 chromatin states ordered by RankDiff boxplot result

out<-"output_directory/"


# Load summary table for epigenetic features
alldf<-read.csv(file = "file_path/E017_summary_StateMean_GenomeMean.csv", header = T)
rownames(alldf)<-alldf[,1]
rownames(alldf)[16]<-"Whole genome"
alldf<-alldf[,2:30]


## Select and subset list of histone marks to keep
colnames(alldf)
Keep<-c(1,2,23,13,17,25,14)
alldf<-alldf[,rev(Keep)]

# Get epigenetic mark names
Names<-colnames(alldf)
EpiNames<-c()
for (i in 1:ncol(alldf)) { 
  Split<-unlist(strsplit(Names[i], "E017."))[2]
  Split<-unlist(strsplit(Split, ".fc.signal.bigwig"))
  EpiNames<-append(EpiNames, Split)
  
}
colnames(alldf)<-EpiNames

# log2 transform fold change values
alldf<-log2(alldf)


## Prepare data for plotting
Newdf<-as.data.frame(alldf)
Newdf$id<-rownames(alldf)
Newdf<-melt(Newdf)
colnames(Newdf)<-c("Chrom_state","Histone_Mark","Log2_FC")
Newdf$Log2_FC[Newdf$value>2]<-2
Newdf$Log2_FC[Newdf$value<(-2)]<-(-2)


# factor levels in order with RankDiff boxplot
Chrom_order<-c("15_Quies","5_TxWk","9_Het","14_ReprPCWk","Whole genome",
               "4_Tx","8_ZNF/Rpts","7_Enh","6_EnhG","1_TssA","3_TxFlnk",
               "13_ReprPC","2_TssAFlnk","10_TssBiv","12_EnhBiv","11_BivFlnk")

Chrom_labels<-c("Quies","TxWk","Het","ReprPCWk","Whole genome",
                "Tx","ZNF/Rpts","Enh","EnhG","TssA","TxFlnk",
                "ReprPC","TssAFlnk","TssBiv","EnhBiv","BivFlnk")

# factor levels
Newdf$Chrom_state<-factor(Newdf$Chrom_state, levels = Chrom_order, labels = Chrom_labels)

# Heatmap plot WITH LABELS
Heatmap_plot<-ggplot(data = Newdf, mapping = aes(x = Chrom_state, y = Histone_Mark,fill = Log2_FC)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient2(name="signal", low = "blue", mid="white", high = "red")+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  ggtitle(label = NULL)+expand_limits(fill=c(-2,2)) +
  theme(plot.title = element_text(size = 10 , face = "bold", hjust = 0.5), 
        axis.text.x = element_text(size = 8, angle = 45, hjust=1), 
        axis.text.y = element_text(size = 8, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
Heatmap_plot
ggsave(filename = "E017_hist_heatmap_chrom_states.RankDiff_order.pdf",
       path = out,plot = Heatmap_plot,device = "pdf",width = 2,height = 2,units = "in")



## END ##
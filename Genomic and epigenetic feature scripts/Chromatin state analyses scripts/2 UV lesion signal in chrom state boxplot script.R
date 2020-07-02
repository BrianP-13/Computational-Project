library(ggplot2)

# Description: This script generates the boxplots showing the distribution of ranked UV lesion signals in IMR90 chromatin states

out<-"output_directory/"


######
## Ranked signal distribution in chromatin states
BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",  "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black","White" )


# 1. Ranked CPD signal
load(file = "file_path/Ranked_CPD_boxplot_stats_for_chrom_states.RData")
CPDdf

CPDdf$BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",
                    "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black","White" )

CPDdf$states<-factor(x = CPDdf$states,levels = c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies","Whole genome"),
                     labels = c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts","Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies", "Whole genome"))

plottitle<-"Ranked CPD signal within 15 Core Chromatin States"
xlab<-"Chromatin States definied in WT IMR90"
ylab<-"CPD (ranked signal)"

# CPD plot with labels
CPD_Plot<-ggplot(CPDdf, aes(x = reorder(states, median))) + 
  geom_boxplot(size=0.25,aes(ymin=min,lower=q25,middle=median,upper=q75,max=max, fill=states),stat="identity") +
  scale_fill_manual(values = BarColors) +
  ggtitle(label = plottitle) + xlab(xlab) + ylab(ylab) + guides(fill=F)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
CPD_Plot
ggsave(filename = "Ranked_CPD_signal_in_chrom_states.pdf",
       path = out,plot = CPD_Plot,device = "pdf",width = 3,height = 2,units = "in")



# 2. Ranked 64PP signal
load(file = "file_path/Ranked_64PP_boxplot_stats_for_chrom_states.RData")
IP64df

IP64df$BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",
                    "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black","White" )
IP64df$states<-factor(x = IP64df$states,levels = c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies","Whole genome"),
                     labels = c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts","Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies", "Whole genome"))
plottitle<-"Ranked 6-4PP signal within 15 Core Chromatin States"

xlab<-"Chromatin States definied in WT IMR90"
ylab<-"6-4PP (ranked signal)"

# 64PP plot with labels
Plot_64<-ggplot(IP64df, aes(x = reorder(states, median))) + 
  geom_boxplot(size=0.25,aes(ymin=min,lower=q25,middle=median,upper=q75,max=max, fill=states),stat="identity") +
  ggtitle(label = plottitle) + xlab(xlab) + ylab(ylab) + 
  scale_fill_manual(values = BarColors) +guides(fill=FALSE)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
Plot_64
ggsave(filename = "Ranked_IP64_signal_in_chrom_states.pdf",
       path = out,plot = Plot_64,device = "pdf",width = 3,height = 2,units = "in")


# 3. Difference
load(file = "file_path/Ranked_64PP_CPD_difference_boxplot_stats_for_chrom_states.RData")
Diffdf

Diffdf$BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",
                     "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black","White" )
Diffdf$states<-factor(x = Diffdf$states,levels = c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies","Whole genome"),
                      labels = c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts","Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies", "Whole genome"))
plottitle<-"Ranked 64PP-CPD within 15 Core Chromatin States"

xlab<-"Chromatin States definied in WT IMR90"
ylab<-"64PP-CPD (ranked signal)"

# Plot with labels
Plot_diff<-ggplot(Diffdf, aes(x = reorder(states, median))) + 
  geom_boxplot(size=0.25,aes(ymin=min,lower=q25,middle=median,upper=q75,max=max, fill=states),stat="identity") +
  scale_fill_manual(values = BarColors) +
  ggtitle(label = plottitle) + xlab(xlab) + ylab(ylab) +guides(fill=FALSE)+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(angle = 45, hjust=1,size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
Plot_diff
ggsave(filename = "Ranked_Diff_in_chrom_states.pdf",
       path = out,plot = Plot_diff,device = "pdf",width = 3,height = 2,units = "in")


## END ##
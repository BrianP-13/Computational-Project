library(ggplot2)
library(data.table)


## Description: Will calculate the sum/cumulative repair score per genomic region. Then I will create a composite plot overlapping
## all repair values for a specific time point per individual chromatin state.


out<-"output_directory/"


# Common labels
core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_labels<-c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts","Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies")
BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",  "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black")
XLab<-"IMR90 - 15 Core Chromatin states"


######
## 64PP

# Load master table
load(file = "file_path/XR64PP_PLUS&MINUS_ALL_TIME_POINTS_core15.RData")
XR64_core15<-XR_table

# Create column with total NER repair score (sum of repair signal at each time point per genomic region)
XR64_core15$Sum_repair<-rowSums(x = XR64_core15[,c(7:11)],dims = 1,na.rm = F)

# Remove NAs
XR64_core15<-na.omit(XR64_core15) 

# Remove sex chromosomes and M chromosome
XR64_core15<-XR64_core15[which(XR64_core15$seqnames!='chrY'),]
XR64_core15<-XR64_core15[which(XR64_core15$seqnames!='chrM'),]
XR64_core15<-XR64_core15[which(XR64_core15$seqnames!='chrX'),]

# Create data.frame for boxplot
XR64_core15<-XR64_core15[,c(6:12)]
colnames(XR64_core15)[1]<-"Chromatin_State"

XR64_core15$BarColors<-0

for (i in 1:15) {
  XR64_core15[which(XR64_core15$Chromatin_State==core15_names[i]),8]<-BarColors[i]
}

# Factor chromatin state labels
XR64_core15$Chromatin_State<-factor(x = XR64_core15$Chromatin_State,levels = core15_names,labels = core15_labels)


## Plot labels
YLab<-"Cumulative 64PP Repair \n (sum of all repair time points)"

# Boxplots with labels
XR64_chrom_plot<-ggplot(data = XR64_core15, aes(x = reorder(Chromatin_State,Sum_repair,FUN = median), y = Sum_repair, fill=Chromatin_State )) +
  geom_boxplot(outlier.shape = NA,show.legend = F,size=0.25) + 
  labs(x= XLab,y= YLab) +
  scale_fill_manual(values = BarColors) + coord_cartesian(ylim = c(25,100))+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust=1,angle = 45))

XR64_chrom_plot
ggsave(filename = "XR64_repair_in_core15_states.pdf",
       path = out,plot = XR64_chrom_plot,device = "pdf",width = 3,height = 2,units = "in")



######
## CPD

# Load master table
load(file = "file_path/XRcpd_PLUS&MINUS_ALL_TIME_POINTS_core15.RData")
XRCPD_core15<-XR_table

# Remove NAs
XRCPD_core15<-na.omit(XRCPD_core15) 

# Create column with total NER repair score (sum of repair signal at each time point per genomic region)
XRCPD_core15$Sum_repair<-rowSums(x = XRCPD_core15[,c(7:12)],dims = 1)

# Remove sex chromosomes and M chromosome
XRCPD_core15<-XRCPD_core15[which(XRCPD_core15$seqnames!='chrY'),]
XRCPD_core15<-XRCPD_core15[which(XRCPD_core15$seqnames!='chrM'),]
XRCPD_core15<-XRCPD_core15[which(XRCPD_core15$seqnames!='chrX'),]

# Create data.frame for boxplot
XRCPD_core15<-XRCPD_core15[,c(6,13)]
colnames(XRCPD_core15)[1]<-"Chromatin_State"

# Add column for color type for each chromatin state
core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
BarColors<-c( "Red","Indian Red","Dark Salmon","Dark Khaki","grey","Green","Dark Green",  "Green Yellow", "Yellow","Medium Aquamarine","Pale Turquoise","Orange Red","Lime Green","Gainsboro","Black")

XRCPD_core15$BarColors<-0

for (i in 1:15) {
  XRCPD_core15[which(XRCPD_core15$Chromatin_State==core15_names[i]),3]<-BarColors[i]
}

XRCPD_core15$Chromatin_State<-factor(x = XRCPD_core15$Chromatin_State,levels = c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies"),labels = c("TssA","TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts","Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"))

# Plot labels
YLab<-"Cumulative CPD Repair \n (sum of all repair time points)"
XLab<-"IMR90 - 15 Core Chromatin states"


# Boxplots with labels
XRCPD_chrom_plot<-ggplot(data = XRCPD_core15, aes(x = reorder(Chromatin_State,Sum_repair,FUN = median), y = Sum_repair, fill=Chromatin_State )) +
  geom_boxplot(outlier.shape = NA,show.legend = F,size=0.25) + 
  labs(x= XLab,y= YLab) +
  scale_fill_manual(values = BarColors) + coord_cartesian(ylim = c(25,120))+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust=1,angle = 45))

XRCPD_chrom_plot
ggsave(filename = "XRcpd_repair_in_core15_states.pdf",
       path = out,plot = XRCPD_chrom_plot,device = "pdf",width = 3,height = 2,units = "in")


## END ##
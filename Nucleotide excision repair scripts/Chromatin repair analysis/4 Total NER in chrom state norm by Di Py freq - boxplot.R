library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
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


# Import dataset with coordinates for 15 core chromatin state 
core15<-file.path('file_path/E017_15_coreMarks_dense.bed.bgz')
df_core15<-df_core15[which(df_core15@seqnames != "chrY"),]
df_core15<-df_core15[which(df_core15@seqnames != "chrX"),]
df_core15<-df_core15[which(df_core15@seqnames != "chrM"),]
mcols(df_core15)$score<-NULL
mcols(df_core15)$itemRgb<-NULL
mcols(df_core15)$thick<-NULL
df_core15<-sortSeqlevels(df_core15)
df_core15<-sort(df_core15)
seqinfo(df_core15)<-hg19
seqlevels(df_core15)<-seqlevelsInUse(df_core15)
df_core15

BSgenome<-BSgenome.Hsapiens.UCSC.hg19
hg19<-seqinfo(BSgenome)

## Calculate diNt freq for all chromatin states
DNAseqs<-getSeq(BSgenome, df_core15)

# Get all di Nt freq for 100kb regions
Freq<-dinucleotideFrequency(DNAseqs, as.prob = F)
Freq<-as.data.frame(Freq)

# Only keep data for dipyrimidines
Freq<-Freq[,c(1,3,6,8,9,11,14,16)] # Include all di-pyrimidine sequences (TT, TC, CT, CC), including on the reverse strand (AA, GA, AG, GG)
Freq$Total_DiPy<-rowSums(Freq[,c(1:8)],dims = 1)

# Add diNt data to chrom state dataset
mcols(df_core15)<-cbind(mcols(df_core15), as.data.frame(Freq))

# Convert back to dataframe object
df_core15<-as.data.frame(df_core15)


# 6-4PP -------------------------------------------------------------------

# Load master table
load(file = "file_path/XR64PP_PLUS&MINUS_ALL_TIME_POINTS_core15.RData")

# Remove sex chromosomes and M chromosome
Total_XR64_Chrom15<-Total_XR64_Chrom15[which(Total_XR64_Chrom15$seqnames!='chrY'),]
Total_XR64_Chrom15<-Total_XR64_Chrom15[which(Total_XR64_Chrom15$seqnames!='chrM'),]
Total_XR64_Chrom15<-Total_XR64_Chrom15[which(Total_XR64_Chrom15$seqnames!='chrX'),]
colnames(Total_XR64_Chrom15)[12]<-"Sum_repair"

## Add total di pyrimidine content data to dataframe
Total_XR64_Chrom15<-cbind(Total_XR64_Chrom15,df_core15[,15])
colnames(Total_XR64_Chrom15)[13]<-"Total_DiPy"

## Normalize total repair by di pyrimidine content
Total_XR64_Chrom15$Norm_total_repair<-Total_XR64_Chrom15$Sum_repair/Total_XR64_Chrom15$Total_DiPy

# Create data.frame for boxplot
XR64_core15<-Total_XR64_Chrom15[,c(4,6,12,13,14)]
colnames(XR64_core15)[2]<-"Chromatin_State"
XR64_core15$BarColors<-0

# Remove regions with zero repair values
XR64_core15<-XR64_core15[which(XR64_core15$Norm_total_repair!=0),]

for (i in 1:15) {
  XR64_core15[which(XR64_core15$Chromatin_State==core15_names[i]),6]<-BarColors[i]
}

# Factor chromatin state labels
XR64_core15$Chromatin_State<-factor(x = XR64_core15$Chromatin_State,levels = core15_names,
                                    labels = core15_labels)

# Plot labels
YLab<-"Normalized cumulative 64PP Repair \n (sum of all repair time points/Total Di-Py counts)"

# Boxplots - Normalized total repair
XR64_chrom_plot_repair<-ggplot(data = XR64_core15, aes(x = reorder(Chromatin_State,Norm_total_repair,FUN = median), y = Norm_total_repair, fill=Chromatin_State )) +
  geom_boxplot(outlier.shape = NA,show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(hjust=1,angle = 45)) + labs(x= XLab,y= YLab) +
  scale_fill_manual(values = BarColors) + coord_cartesian(ylim = c(0,0.2))
XR64_chrom_plot_repair

XR64_chrom_plot_name<-paste0(out,"IP64/","XR64_repair_in_core15_states_DiPy_normalized.jpeg")
jpeg(filename = XR64_chrom_plot_name,width = 4,height = 4,units = "in",res = 300)
print(XR64_chrom_plot_repair)
dev.off()



# CPD ---------------------------------------------------------------------

# Load master table
load(file = "file_path/Total_XRCPD_strand_mean_repair_core15.RData")

# Remove sex chromosomes and M chromosome
Total_XRCPD_Chrom15<-Total_XRCPD_Chrom15[which(Total_XRCPD_Chrom15$seqnames!='chrY'),]
Total_XRCPD_Chrom15<-Total_XRCPD_Chrom15[which(Total_XRCPD_Chrom15$seqnames!='chrM'),]
Total_XRCPD_Chrom15<-Total_XRCPD_Chrom15[which(Total_XRCPD_Chrom15$seqnames!='chrX'),]
colnames(Total_XRCPD_Chrom15)[13]<-"Sum_repair"

## Add total di pyrimidine content data to dataframe
Total_XRCPD_Chrom15<-cbind(Total_XRCPD_Chrom15,df_core15[,15])
colnames(Total_XRCPD_Chrom15)[14]<-"Total_DiPy"

## Normalize total repair by di pyrimidine content
Total_XRCPD_Chrom15$Norm_total_repair<-Total_XRCPD_Chrom15$Sum_repair/Total_XRCPD_Chrom15$Total_DiPy

# Create data.frame for boxplot
XRcpd_core15<-Total_XRCPD_Chrom15[,c(4,6,13,14,15)]
colnames(XRcpd_core15)[2]<-"Chromatin_State"
XRcpd_core15$BarColors<-0

# Remove regions with zero repair values
XRcpd_core15<-XRcpd_core15[which(XRcpd_core15$Norm_total_repair!=0),]

for (i in 1:15) {
  XRcpd_core15[which(XRcpd_core15$Chromatin_State==core15_names[i]),6]<-BarColors[i]
}

# Factor chromatin state labels
XRcpd_core15$Chromatin_State<-factor(x = XRcpd_core15$Chromatin_State,levels = core15_names,
                                    labels = core15_labels)
# Plot labels
YLab<-"Normalized cumulative CPD Repair \n (sum of all repair time points/Total Di-Py counts)"

# Boxplots
XRcpd_chrom_plot<-ggplot(data = XRcpd_core15, aes(x = reorder(Chromatin_State,Norm_total_repair,FUN = median), y = Norm_total_repair, fill=Chromatin_State )) +
  geom_boxplot(outlier.shape = NA,show.legend = F) + theme_bw() + 
  theme(axis.text.x = element_text(hjust=1,angle = 45)) + labs(x= XLab,y= YLab) +
  scale_fill_manual(values = BarColors) + coord_cartesian(ylim = c(0,0.25))
XRcpd_chrom_plot

XRcpd_chrom_plot_name<-paste0(out,"CPD/","XRcpd_repair_in_core15_states_DiPy_normalized.jpeg")
jpeg(filename = XRcpd_chrom_plot_name,width = 4,height = 4,units = "in",res = 300)
print(XRcpd_chrom_plot)
dev.off()


### END ###
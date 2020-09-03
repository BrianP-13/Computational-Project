library(GenomicRanges)
library(reshape2)
library(IDPmisc)
library(ggplot2)
library(biomaRt)


out<-"output_directory/"


# ggplot themes
THEMES_1<-theme(plot.title = element_text(size = 8 , face = "bold", hjust = 0.5), 
                axis.text.x = element_text(size = 6, angle = 45, hjust=1), 
                axis.text.y = element_text(size = 6, hjust=1),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank())
THEMES_2<-theme(line = element_line(size = 0.15), 
                panel.border = element_rect(fill = NA,size = 0.25),
                panel.grid.minor = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank())

#
# Load melanoma cancer driver list with mutation frequency info
load(file = 'file_path/MELA_driver_list_with_mut_freq_info.RData')
MELA_cancer_drivers


# Load dataset with all UV lesion info and mutation rate per gene
load('file_path/Prot_coding_genes_w_repair_mutR_UVratio_info.RData')
UV_dmg_repair_mutR<-subset(UV_dmg_repair_mutR,seqnames!="chrX")

### SETTING UP LIMITS FOR ALL VARIABLE (REPAIR, MUTATION RATE)
# Set upper limit of repair values

# CPD damage
CPD_dmg_upper_limit<-boxplot.stats(UV_dmg_repair_mutR$Log2_CPD)$stats[5]
CPD_dmg_mid_point<-boxplot.stats(UV_dmg_repair_mutR$Log2_CPD)$stats[3]
CPD_dmg_lower_limit<-boxplot.stats(UV_dmg_repair_mutR$Log2_CPD)$stats[1]

# CPD repair
XRcpd_upper_limit<-boxplot.stats(UV_dmg_repair_mutR[which(UV_dmg_repair_mutR$Total_XRCPD_strand_mean_repair!=0),13])$stats[5]
XRcpd_mid_point<-boxplot.stats(UV_dmg_repair_mutR[which(UV_dmg_repair_mutR$Total_XRCPD_strand_mean_repair!=0),13])$stats[3]

# Set upper and lower limit of mutation rate values
MutR_upper_limit<-boxplot.stats(UV_dmg_repair_mutR$Log10_mut_rate)$stats[5]
MutR_med<-boxplot.stats(UV_dmg_repair_mutR$Log10_mut_rate)$stats[3]
MutR_lower_limit<-min(UV_dmg_repair_mutR[which(UV_dmg_repair_mutR$Log10_mut_rate!="-Inf"),15])


######
# Top CPD ratio
Top_10_CPDratio_val<-quantile(UV_dmg_repair_mutR$Log2_UVcpd_ratio, 0.9, na.rm=T)
CPD_ratio_T10<-UV_dmg_repair_mutR[which(UV_dmg_repair_mutR$Log2_UVcpd_ratio >=Top_10_CPDratio_val ),]
CPD_ratio_T10$Region<-"Top_10_CPD"

CM_drivers_Tcpd<-CPD_ratio_T10[which(CPD_ratio_T10$gene_id %in% MELA_cancer_drivers$gene_id),]


# Change values in dataset with new limits
CM_drivers_Tcpd$Log2_CPD[CM_drivers_Tcpd$Log2_CPD>CPD_dmg_upper_limit]<-CPD_dmg_upper_limit
CM_drivers_Tcpd$Log2_CPD[CM_drivers_Tcpd$Log2_CPD<CPD_dmg_lower_limit]<-CPD_dmg_lower_limit

CM_drivers_Tcpd$Total_XRCPD_strand_mean_repair[CM_drivers_Tcpd$Total_XRCPD_strand_mean_repair>XRcpd_upper_limit]<-XRcpd_upper_limit

CM_drivers_Tcpd$Log10_mut_rate[CM_drivers_Tcpd$Log10_mut_rate>MutR_upper_limit]<-MutR_upper_limit
CM_drivers_Tcpd$Log10_mut_rate[CM_drivers_Tcpd$Log10_mut_rate<MutR_lower_limit]<-MutR_lower_limit


######
## UV damage heatmaps

# Heatmap of CPD damage
TCPD_CPDdamage_heatmap<-ggplot(data = CM_drivers_Tcpd, mapping = aes(x = Region, y = hgnc_symbol,fill = Log2_CPD)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "red", 
                      limits=c(CPD_dmg_lower_limit,CPD_dmg_upper_limit),
                      breaks=c(CPD_dmg_lower_limit,CPD_dmg_mid_point,CPD_dmg_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
TCPD_CPDdamage_heatmap
ggsave(filename = "TCPD_ratio_driver_genes_CPDdmg_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = TCPD_CPDdamage_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


######
## UV repair heatmaps

# Heatmap of CPD repair
TCPD_CPDrepair_heatmap<-ggplot(data = CM_drivers_Tcpd, mapping = aes(x = Region, y = hgnc_symbol,fill = Total_XRCPD_strand_mean_repair)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "green", 
                      limits=c(0,XRcpd_upper_limit), breaks=c(0,XRcpd_mid_point,XRcpd_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  ggtitle(label = NULL)+
  expand_limits(fill=c(0,XRcpd_upper_limit))+THEMES_1
TCPD_CPDrepair_heatmap
ggsave(filename = "TCPD_ratio_driver_genes_XRCPD_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = TCPD_CPDrepair_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


######
## mutation rate heatmaps

# Heatmap of mutation rate
TCPD_MutR_heatmap<-ggplot(data = CM_drivers_Tcpd, mapping = aes(x = Region, y = hgnc_symbol,fill = Log10_mut_rate)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "yellow", 
                      limits=c(MutR_lower_limit,MutR_upper_limit),breaks=c(MutR_lower_limit,MutR_med,MutR_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
TCPD_MutR_heatmap
ggsave(filename = "TCPD_ratio_driver_genes_MutR_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = TCPD_MutR_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


### END ###
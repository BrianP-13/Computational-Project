library(GenomicRanges)
library(reshape2)
library(IDPmisc)
library(ggplot2)



out<-"output_directory/"


# Overlap function:
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


# Load enhancer gene interactions
load(file = 'file_path/IMR90_Enh_prot_coding_gene_interactions.RData')
IMR90_Enh_prot_genes<-as.data.frame(IMR90_Enh_prot_genes)
IMR90_Enh_prot_genes<-subset(IMR90_Enh_prot_genes,seqnames!="chrX")
length(unique(IMR90_Enh_prot_genes$Gene_name)) # total of 11,294 genes with enhancer pair

# Load enhancer dataset with lesion, repair, UV ratio and mutation rate info
load('file_path/Enhancers_w_repair_mutR_UVratio_info.RData')


## Need to set up limits for all values

### SETTING UP LIMITS FOR ALL VARIABLE (REPAIR, MUTATION RATE)
# Set upper limit of repair values

# 64PP damage
UV64_data<-Enh_UVdmg_repair_mutR[,c(1:7,9:11,13,15,16)]
UV64_data<-subset(UV64_data,UV64_ratio>=0) # Only include enhancers that have lesion susc and repair values

IP64_dmg_upper_limit<-boxplot.stats(UV64_data$Log2_64PP)$stats[5]
IP64_dmg_mid_point<-boxplot.stats(UV64_data$Log2_64PP)$stats[3]
IP64_dmg_lower_limit<-boxplot.stats(UV64_data$Log2_64PP)$stats[1]

# 64PP repair
XR64_upper_limit<-boxplot.stats(UV64_data[which(UV64_data$Total_XR64_strand_mean_repair!=0),7])$stats[5]
XR64_mid_point<-boxplot.stats(UV64_data[which(UV64_data$Total_XR64_strand_mean_repair!=0),7])$stats[3]

#
# CPD damage
UVcpd_data<-Enh_UVdmg_repair_mutR[,c(1:6,8:10,12,14,17,18)]
UVcpd_data<-subset(UVcpd_data,UVcpd_ratio>=0) # Only include enhancers that have lesion susc and repair values

CPD_dmg_upper_limit<-boxplot.stats(UVcpd_data$Log2_CPD)$stats[5]
CPD_dmg_mid_point<-boxplot.stats(UVcpd_data$Log2_CPD)$stats[3]
CPD_dmg_lower_limit<-boxplot.stats(UVcpd_data$Log2_CPD)$stats[1]

# CPD repair
XRcpd_upper_limit<-boxplot.stats(UVcpd_data[which(UVcpd_data$Total_XRCPD_strand_mean_repair!=0),7])$stats[5]
XRcpd_mid_point<-boxplot.stats(UVcpd_data[which(UVcpd_data$Total_XRCPD_strand_mean_repair!=0),7])$stats[3]

# Set upper and lower limit of mutation rate values
MutR_upper_limit<-boxplot.stats(UVcpd_data$Log10_mut_rate)$stats[5]
MutR_lower_limit<-boxplot.stats(UVcpd_data$Log10_mut_rate)$stats[1]


######
# Top CPD ratio
Top_10_CPDratio_val<-quantile(Enh_UVdmg_repair_mutR$Log2_UVcpd_ratio, 0.9, na.rm=T)
CPD_ratio_T10<-Enh_UVdmg_repair_mutR[which(Enh_UVdmg_repair_mutR$Log2_UVcpd_ratio >=Top_10_CPDratio_val),]

# Get list of enhancer-genes that overlap with enhancers with top UV ratio scores
T10_CPD_Enh<-Ovrlp_fun(Query = makeGRangesFromDataFrame(IMR90_Enh_prot_genes,keep.extra.columns = T), Subject = makeGRangesFromDataFrame(CPD_ratio_T10,keep.extra.columns = T), Overlap = 50, Region_name = "Top 10% CPD ratio scrore")
T10_CPD_Enh<-as.data.frame(T10_CPD_Enh)

# Subset list of enhancer-genes to only include cancer driver genes
CM_drivers_Tcpd<-T10_CPD_Enh[which(T10_CPD_Enh$Gene_ID %in% MELA_cancer_drivers$gene_id),]

# Create list of driver genes in selected group
Tcpd_unique_drivers<-unique(CM_drivers_Tcpd$Gene_name)

# Set up new table to save averages
Tcpd_Enh_driver_df<-as.data.frame( matrix(ncol = 15,nrow = 0))


for (i in 1:length(Tcpd_unique_drivers)) {
  Tcpd_unique_drivers[i]
  df<-CM_drivers_Tcpd[which(CM_drivers_Tcpd$Gene_name==Tcpd_unique_drivers[i]),]
  
  enhs<-subsetByOverlaps(x = makeGRangesFromDataFrame(CPD_ratio_T10,keep.extra.columns = T),
                         ranges = makeGRangesFromDataFrame(df,keep.extra.columns = T))
  
  enhs_df<-as.data.frame(enhs)
  enhs_df<-enhs_df[,c(1:6,12,8,9,17,18)]
  
  enhs_df$Enh_target_gene_name<-Tcpd_unique_drivers[i]
  enhs_df$'Mean XRCPD signal(if mult enhancers)'<-NA
  enhs_df$'Mean CPD FC signal(if mult enhancers)'<-NA
  enhs_df$'Mean C>T rate (if mult enhancers)'<-NA
  
  enhs_df[1,13]<-mean(enhs_df$Total_XRCPD_strand_mean_repair)
  enhs_df[1,14]<-mean(enhs_df$Mean_CPD_FC)
  enhs_df[1,15]<-mean(enhs_df$Total_C_to_T_rate)
  
  colnames(Tcpd_Enh_driver_df)<-colnames(enhs_df)
  Tcpd_Enh_driver_df<-rbind(Tcpd_Enh_driver_df,enhs_df)
  
}
Tcpd_Enh_driver_df

write.csv(Tcpd_Enh_driver_df, file = paste0(out,"Driver genes reg by Enh with top CPD ratio scores.csv"))


# Make new dataframe for heatmap plotting
Tcpd_Enh_driver_df2<-Tcpd_Enh_driver_df[,c(12:15)]
rownames(Tcpd_Enh_driver_df2)<-NULL
Tcpd_Enh_driver_df2<-na.omit(Tcpd_Enh_driver_df2)

colnames(Tcpd_Enh_driver_df2)[2]<-"Total_XRCPD_strand_mean_repair"
Tcpd_Enh_driver_df2$Log2_CPD<-log2(Tcpd_Enh_driver_df2$`Mean CPD FC signal(if mult enhancers)`)
Tcpd_Enh_driver_df2$Log10_mut_rate<-log10(Tcpd_Enh_driver_df2$`Mean C>T rate (if mult enhancers)`)


## Change values in dataset with new limits
# UV susceptibility limits
Tcpd_Enh_driver_df2$Log2_CPD[Tcpd_Enh_driver_df2$Log2_CPD>CPD_dmg_upper_limit]<-CPD_dmg_upper_limit
Tcpd_Enh_driver_df2$Log2_CPD[Tcpd_Enh_driver_df2$Log2_CPD<CPD_dmg_lower_limit]<-CPD_dmg_lower_limit

# UV repair limits
Tcpd_Enh_driver_df2$Total_XRCPD_strand_mean_repair[Tcpd_Enh_driver_df2$Total_XRCPD_strand_mean_repair>XRcpd_upper_limit]<-XRcpd_upper_limit

# Mutation rate limits
MutR_upper_limit<-boxplot.stats(UVcpd_data$Log10_mut_rate)$stats[5]
MutR_med<-boxplot.stats(UVcpd_data$Log10_mut_rate)$stats[3]
MutR_lower_limit<-boxplot.stats(UVcpd_data$Log10_mut_rate)$stats[1]
MutR_lower_limit<-min(UVcpd_data$Log10_mut_rate)

Tcpd_Enh_driver_df2$Log10_mut_rate[Tcpd_Enh_driver_df2$Log10_mut_rate>MutR_upper_limit]<-MutR_upper_limit
Tcpd_Enh_driver_df2$Log10_mut_rate[Tcpd_Enh_driver_df2$Log10_mut_rate<MutR_lower_limit]<-MutR_lower_limit


## Heatmap plots
Tcpd_Enh_driver_df2$Region<-"Top 10% CPD\nratio scrore"

######
## UV damage heatmaps

# Heatmap of 64PP damage
Tcpd_CPDdamage_heatmap<-ggplot(data = Tcpd_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Log2_CPD)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "red", 
                      limits=c(CPD_dmg_lower_limit,CPD_dmg_upper_limit),
                      breaks=c(CPD_dmg_lower_limit,CPD_dmg_mid_point,CPD_dmg_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
Tcpd_CPDdamage_heatmap
ggsave(filename = "Enh_Tcpd_ratio_driver_genes_CPDdmg_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = Tcpd_CPDdamage_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


######
## UV repair heatmaps

# Heatmap of 64PP repair
Tcpd_CPDrepair_heatmap<-ggplot(data = Tcpd_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Total_XRCPD_strand_mean_repair)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "green", 
                      limits=c(0,XRcpd_upper_limit), breaks=c(0,XRcpd_mid_point,XRcpd_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
Tcpd_CPDrepair_heatmap
ggsave(filename = "Enh_Tcpd_ratio_driver_genes_XRCPD_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = Tcpd_CPDrepair_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")



######
## Heatmap of mutation rate

Tcpd_MutR_heatmap<-ggplot(data = Tcpd_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Log10_mut_rate)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "yellow", 
                      limits=c(MutR_lower_limit,MutR_upper_limit),breaks=c(MutR_lower_limit,MutR_med,MutR_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
Tcpd_MutR_heatmap
ggsave(filename = "Enh_Tcpd_ratio_driver_genes_MutR_heatmap.pdf",
       path = paste0(out, "CPD/"),plot = Tcpd_MutR_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")




### END ###
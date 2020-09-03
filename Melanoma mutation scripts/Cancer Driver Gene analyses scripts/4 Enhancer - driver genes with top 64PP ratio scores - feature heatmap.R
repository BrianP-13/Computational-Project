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
MELA_cancer_drivers

# Load enhancer gene interactions
load(file = 'file_path/IMR90_Enh_prot_coding_gene_interactions.RData')
IMR90_Enh_prot_genes<-as.data.frame(IMR90_Enh_prot_genes)
IMR90_Enh_prot_genes<-subset(IMR90_Enh_prot_genes,seqnames!="chrX")
length(unique(IMR90_Enh_prot_genes$Gene_name)) # total of 11,294 genes with enhancer pair


# Load enhancer dataset with lesion, repair, UV ratio and mutation rate info
load('file_path/Enhancers_w_repair_mutR_UVratio_info.RData')


### SETTING UP LIMITS FOR ALL VARIABLE (REPAIR, MUTATION RATE)
# Set upper limit of repair values

# 64PP damage
UV64_data<-Enh_UVdmg_repair_mutR[,c(1:7,9:11,13,15,16)]
UV64_data<-subset(UV64_data,UV64_ratio>=0) # Only include enhancers that have lesion susc and repair values

IP64_dmg_upper_limit<-boxplot.stats(UV64_data$Log2_64PP)$stats[5]
IP64_dmg_mid_point<-boxplot.stats(UV64_data$Log2_64PP)$stats[3]
IP64_dmg_lower_limit<-boxplot.stats(UV64_data$Log2_64PP)$stats[1]

# 64PP repair
boxplot.stats(UV64_data[which(UV64_data$Total_XR64_strand_mean_repair!=0),7])$stats
XR64_upper_limit<-boxplot.stats(UV64_data[which(UV64_data$Total_XR64_strand_mean_repair!=0),7])$stats[5]
XR64_mid_point<-boxplot.stats(UV64_data[which(UV64_data$Total_XR64_strand_mean_repair!=0),7])$stats[3]

# Set upper and lower limit of mutation rate values
MutR_upper_limit<-boxplot.stats(UV64_data$Log10_mut_rate)$stats[5]
MutR_med<-boxplot.stats(UV64_data$Log10_mut_rate)$stats[3]
MutR_lower_limit<-min(UV64_data$Log10_mut_rate)


# Histograms
head(UV64_data)
ggplot(data = UV64_data, aes(x = Total_XR64_strand_mean_repair))+geom_histogram()+
  geom_vline(xintercept = XR64_upper_limit,color="red")+
  geom_vline(xintercept = XR64_mid_point,color="green")


######
# Top 64PP ratio
Top_10_64ratio_val<-quantile(Enh_UVdmg_repair_mutR$Log2_UV64_ratio, 0.9, na.rm=T)
IP64_ratio_T10<-Enh_UVdmg_repair_mutR[which(Enh_UVdmg_repair_mutR$Log2_UV64_ratio >=Top_10_64ratio_val),]

# Get list of enhancer-genes that overlap with enhancers with top UV ratio scores
T10_64_Enh<-Ovrlp_fun(Query = makeGRangesFromDataFrame(IMR90_Enh_prot_genes,keep.extra.columns = T), Subject = makeGRangesFromDataFrame(IP64_ratio_T10,keep.extra.columns = T), Overlap = 50, Region_name = "Top 10% 64PP ratio scrore")
T10_64_Enh<-as.data.frame(T10_64_Enh)

# Subset list of enhancer-genes to only include cancer driver genes
CM_drivers_T64<-T10_64_Enh[which(T10_64_Enh$Gene_ID %in% MELA_cancer_drivers$gene_id),]

# Create list of driver genes in selected group
T64_unique_drivers<-unique(CM_drivers_T64$Gene_name)

# Set up new table to save averages
T64_Enh_driver_df<-as.data.frame( matrix(ncol = 15,nrow = 0))

for (i in 1:length(T64_unique_drivers)) {
  T64_unique_drivers[i]
  df<-CM_drivers_T64[which(CM_drivers_T64$Gene_name==T64_unique_drivers[i]),]
  
  enhs<-subsetByOverlaps(x = makeGRangesFromDataFrame(IP64_ratio_T10,keep.extra.columns = T),
                         ranges = makeGRangesFromDataFrame(df,keep.extra.columns = T))
  
  enhs_df<-as.data.frame(enhs)
  enhs_df<-enhs_df[,c(1:6,11,7,9,15,16)]
  
  enhs_df$Enh_target_gene_name<-T64_unique_drivers[i]
  enhs_df$'Mean XR64 signal(if mult enhancers)'<-NA
  enhs_df$'Mean 64PP FC signal(if mult enhancers)'<-NA
  enhs_df$'Mean C>T rate (if mult enhancers)'<-NA
  
  enhs_df[1,13]<-mean(enhs_df$Total_XR64_strand_mean_repair)
  enhs_df[1,14]<-mean(enhs_df$Mean_64PP_FC)
  enhs_df[1,15]<-mean(enhs_df$Total_C_to_T_rate)
  
  colnames(T64_Enh_driver_df)<-colnames(enhs_df)
  T64_Enh_driver_df<-rbind(T64_Enh_driver_df,enhs_df)
  
}
T64_Enh_driver_df

write.csv(T64_Enh_driver_df, file = paste0(out,"Driver genes reg by Enh with top 64PP ratio scores.csv"))


# Make new dataframe for heatmap plotting
T64_Enh_driver_df2<-T64_Enh_driver_df[,c(12:15)]
rownames(T64_Enh_driver_df2)<-NULL
T64_Enh_driver_df2<-na.omit(T64_Enh_driver_df2)

colnames(T64_Enh_driver_df2)[2]<-"Total_XR64_strand_mean_repair"
T64_Enh_driver_df2$Log2_64PP<-log2(T64_Enh_driver_df2$`Mean 64PP FC signal(if mult enhancers)`)
T64_Enh_driver_df2$Log10_mut_rate<-log10(T64_Enh_driver_df2$`Mean C>T rate (if mult enhancers)`)

## Change values in dataset with new limits
# UV susceptibility limits
T64_Enh_driver_df2$Log2_64PP[T64_Enh_driver_df2$Log2_64PP>IP64_dmg_upper_limit]<-IP64_dmg_upper_limit
T64_Enh_driver_df2$Log2_64PP[T64_Enh_driver_df2$Log2_64PP<IP64_dmg_lower_limit]<-IP64_dmg_lower_limit

# UV repair limits
T64_Enh_driver_df2$Total_XR64_strand_mean_repair[T64_Enh_driver_df2$Total_XR64_strand_mean_repair>XR64_upper_limit]<-XR64_upper_limit

# Mutation rate limits
T64_Enh_driver_df2$Log10_mut_rate[T64_Enh_driver_df2$Log10_mut_rate>MutR_upper_limit]<-MutR_upper_limit
T64_Enh_driver_df2$Log10_mut_rate[T64_Enh_driver_df2$Log10_mut_rate<MutR_lower_limit]<-MutR_lower_limit
head(T64_Enh_driver_df2)


## Heatmap plots
T64_Enh_driver_df2$Region<-"Top 10% 64PP\nratio scrore"

######
## UV damage heatmaps

# Heatmap of 64PP damage
T64_64damage_heatmap<-ggplot(data = T64_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Log2_64PP)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "red", 
                      limits=c(IP64_dmg_lower_limit,IP64_dmg_upper_limit),
                      breaks=c(IP64_dmg_lower_limit,IP64_dmg_mid_point,IP64_dmg_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
T64_64damage_heatmap
ggsave(filename = "Enh_T64_ratio_driver_genes_64dmg_heatmap.pdf",
       path = paste0(out, "IP64/"),plot = T64_64damage_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


######
## UV repair heatmaps

# Heatmap of 64PP repair
T64_64repair_heatmap<-ggplot(data = T64_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Total_XR64_strand_mean_repair)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "green", 
                      limits=c(0,XR64_upper_limit), breaks=c(0,XR64_mid_point,XR64_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
T64_64repair_heatmap
ggsave(filename = "Enh_T64_ratio_driver_genes_XR64_heatmap.pdf",
       path = paste0(out, "IP64/"),plot = T64_64repair_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")


######
## Heatmap of mutation rate
T64_MutR_heatmap<-ggplot(data = T64_Enh_driver_df2, mapping = aes(x = Region, y = Enh_target_gene_name,fill = Log10_mut_rate)) +
  geom_tile() + theme_bw()+
  scale_fill_gradient(name="signal", low = "black", high = "yellow", 
                      limits=c(MutR_lower_limit,MutR_upper_limit),breaks=c(MutR_lower_limit,MutR_med,MutR_upper_limit))+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
  #expand_limits(fill=c(0,XR64_upper_limit)) +
  ggtitle(label = NULL)+THEMES_1
T64_MutR_heatmap
ggsave(filename = "Enh_T64_ratio_driver_genes_MutR_heatmap.pdf",
       path = paste0(out, "IP64/"),plot = T64_MutR_heatmap,device = "pdf",width = 1.5,height = 4,units = "in")




### END ###
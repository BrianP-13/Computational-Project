library(ggplot2)


## Description: This script is for generating a heatmap to show enrichment results of genic features and enhnacers in 15 chromatin states


out<-"output_directory/"


Rep_step<-function(df) {
  df<-df[c(3,6),]
  df<-as.data.frame(t(df))
  
  df<-log10(df)
  df[df <= 0]<-0
  
  df$Diff<- df[,1]-df[,2]
  df$feature<-as.character(rownames(df))
  return(df)
}

# function with yellow/blue combinations
Gplot<-function(Data, gtitle, aes_x, aes_y,aes_fill,Leg) {
  Plot<-ggplot(data = Data, mapping = aes_string(x = aes_x, y = aes_y,fill=aes_fill)) +
    geom_tile(color="darkgrey",size=0.15, show.legend = Leg) + 
    scale_fill_gradient2(name="p-value", low = "dodgerblue", mid="black", high = "yellow",labels=expression(10^6,10^3, 0,10^3,10^6))+
    ggtitle(label = gtitle) +
    theme(plot.title = element_text(size = 10 , face = "bold", hjust = 0.5), 
          axis.text.x = element_text(size = 8, angle = 45, hjust=1), 
          axis.text.y = element_text(size = 8, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  return(Plot)
  
}


# Load enrichment results
load(file = "file_path/Genic.EnhAltas_15ChromState_Enrchment.RData")
str(Enrich_list)

core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_names2<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF_Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")

# Annotations
Annots<-c("Super_enhancer","Enhancer-EnhancerAtlas","hg19_genes_1to5kb","hg19_genes_promoters", "hg19_genes_5UTRs", "hg19_genes_exons", "hg19_genes_introns", "hg19_genes_3UTRs","hg19_genes_intergenic")
Annots2<-c("super \n enhancers","enhancers", "1 to 5kb","promoters","5'UTRs", "exons", "introns", "3'UTRs", "intergenic")


#####
## Order of chromatin states to match distribution of RankDiff

# Format and extract enrichment results from list:
GenicF_df<-as.data.frame(matrix(ncol = 3,nrow = 0))

for (i in 1:15) {
  df<- Enrich_list[[i]]
  df<-Rep_step(df = df)
  df$Chrom_state<-core15_names2[i]
  
  GenicF_df<-rbind(GenicF_df,df,make.row.names=F)
}

# Add extra Whole Genome sample 
WG_df<- Enrich_list[[1]]
WG_df<-Rep_step(df = WG_df)
WG_df$Diff<-0
WG_df$Chrom_state<-"Genome median"
GenicF_df<-rbind(GenicF_df,WG_df,make.row.names=F)

GenicF_df$Diff[GenicF_df$Diff >= 6]= 6
GenicF_df$Diff[GenicF_df$Diff <= -6]= -6

# Add whole genome dataset and change order of chrom states
core15_names3<-append(core15_names2,"Genome median")
RankDiff_order<-core15_names3[c(15,5,9,14,16,4,8,7,6,1,3,13,2,10,12,11)]

# Add factor levels
GenicF_df$Chrom_state<-factor(GenicF_df$Chrom_state, levels = RankDiff_order,
                                labels = c("Quies","TxWk","Het","ReprPCWk","Genome median","Tx","ZNF/Rpts","Enh","EnhG",
                                           "TssA","TxFlnk","ReprPC","TssAFlnk","TssBiv","EnhBiv","BivFlnk"))
GenicF_df$feature<-factor(GenicF_df$feature,levels = rev(Annots),
                          labels =rev(Annots2))

# Heatmap with labels and legend
GenicF_RankDiff_plot<-Gplot(Data = GenicF_df,gtitle = NULL,aes_x = 'Chrom_state',aes_y = 'feature',aes_fill = 'Diff',Leg = T)+
  scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))
GenicF_RankDiff_plot
ggsave(filename = "GenicF_EnhAtlas_Enrch_15core.RANKDIFF.pdf",
       path = out,plot = GenicF_RankDiff_plot,device = "pdf",width = 3,height = 2,units = "in")



## END ##
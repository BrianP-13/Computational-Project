library(ggplot2)
library(IDPmisc)


## Description: Plot with all genic and enhancer elements from EnhancerAtlas database from all regions of interest.

out<-"output_directory/"

# Function to process the enrichment table for ggplot mapping
Rep_step<-function(df) {
  df<-df[c(3,6),]
  df<-as.data.frame(t(df))
  
  df<-log10(df)
  df[df <= 0]<-0
  
  df$Diff<- df[,1]-df[,2]
  df$feature<-as.character(rownames(df))
  df$Color<- ifelse(df$Diff < 0, "negative","positive")
  return(df)
}

load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

# Annotations
SupEnh<-"Super\nenhancers"
Enh<-"EnhAtlas\nenhancers"
genicF<-c("promoters","1 to 5kb", "5'UTRs", "exons", "introns", "3'UTRs", "intergenic")
All_features<-c("Super\nenhancers","EnhAtlas\nenhancers", "1 to 5kb", "promoters", "5'UTRs", "exons", "introns", "3'UTRs", "intergenic")

# Common plot labels
xLab<-"Enhancer and genic features"
yLab = bquote(atop(bold("Depletion <--> Enrichment"), "-log10(pval)" ))


## File names in list:
#  Enh_enrch
Enh_Top_Names<-c("IP64_T10_EnhAtlas_EnhElmts_enrich_100kb.RData","CPD_T10_EnhAtlas_EnhElmts_enrich_100kb.RData","RankDiff_T10_EnhAtlas_EnhElmts_enrich_100kb.RData")
Enh_Bottom_Names<-c("IP64_B10_EnhAtlas_EnhElmts_enrich_100kb.RData","CPD_B10_EnhAtlas_EnhElmts_enrich_100kb.RData","RankDiff_B10_EnhAtlas_EnhElmts_enrich_100kb.RData")

#  Super_Enh_enrch
Sup_Enh_Top_Names<-c("IP64_T10_SUPER_Enh_enrich_100kb.RData","CPD_T10_SUPER_Enh_enrich_100kb.RData","RankDiff_T10_SUPER_Enh_enrich_100kb.RData")
Sup_Enh_Bottom_Names<-c("IP64_B10_SUPER_Enh_enrich_100kb.RData","CPD_B10_SUPER_Enh_enrich_100kb.RData","RankDiff_B10_SUPER_Enh_enrich_100kb.RData")

# Genic_enrch
Genic_Names_Top<-c("IP64_T10_GenicF_enrich_100kb.RData","CPD_T10_GenicF_enrich_100kb.RData","RankDiff_T10_GenicF_enrich_100kb.RData")
Genic_Names_Bottom<-c("IP64_B10_GenicF_enrich_100kb.RData","CPD_B10_GenicF_enrich_100kb.RData","RankDiff_B10_GenicF_enrich_100kb.RData")

# Sample/group names
Sample_Names<-c("64PP","CPD","RankDiff")
Sample_Names_Top<-c("Top_64PP","Top_CPD","Top_RankDiff")
Sample_Names_Bottom<-c("Bottom_64PP","Bottom_CPD","Bottom_RankDiff")

path<-"file_path/"


######
## Top regions:
All_df<-as.data.frame(matrix(nrow = 0,ncol = 6)  )

for (i in 1:3){
  Enh_df<-load_obj(f = paste0(path,Enh_Top_Names[i]))
  SupEnh_df<-load_obj(f = paste0(path,Sup_Enh_Top_Names[i]))
  Genic_df<-load_obj(f = paste0(path,Genic_Names_Top[i]))
  
  colnames(Enh_df)<-Enh
  colnames(SupEnh_df)<-SupEnh
  colnames(Genic_df)<-genicF
  Genic_df<-Genic_df[,c(2,1,3:7)]
  
  G_E_df<-cbind(SupEnh_df,Enh_df,Genic_df)
  G_E_df<-Rep_step(G_E_df)
  G_E_df$feature<-factor(x = G_E_df$feature,levels = All_features)
  G_E_df$Group<-Sample_Names[i]
  
  All_df<-rbind(All_df, G_E_df)
  colnames(All_df)<-colnames(G_E_df)
  
}
All_df$Group<-factor(x = All_df$Group,levels = Sample_Names)
All_df$Diff[All_df$Diff >= 6]=6
All_df$Diff[All_df$Diff <= -6]=-6

# TOP regions
Plot_Top<-ggplot(data = All_df, aes(x=feature, y=Diff), color=Group) + 
  geom_bar(stat = "identity", aes(fill=Group), show.legend = F,width = 0.75,position = "dodge") +theme_bw() +
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),
                     labels = expression(10^6,10^3, 0,10^3,10^6),
                     minor_breaks = NULL) +
  scale_x_discrete(labels= All_features)+
  theme(plot.title = element_text(size = 10, face = "bold" ,hjust = 0.5),
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 45,hjust = 1),axis.ticks.y.left = element_line(size = 0.25))+
  labs(x=xLab,y=yLab, title = NULL, parse=T)
Plot_Top
ggsave(filename = "ALL_T10.EnhAtlas_EnhE.GenicF.enrich.100kb.pdf",
       path = out,plot = Plot_Top,device = "pdf",width = 2,height = 1,units = "in")

# TOP regions- plot rotated 90 degrees, no labels:
Plot_Top2<-ggplot(data = All_df, aes(x=feature, y=Diff), color=Group) + 
  geom_vline(xintercept = 1.5,color="light gray",size = 0.15)+geom_vline(xintercept = 2.5,color="light gray",size = 0.15)+geom_vline(xintercept = 3.5,color="light gray",size = 0.15)+geom_vline(xintercept = 4.5,color="light gray",size = 0.15)+geom_vline(xintercept = 5.5,color="light gray",size = 0.15)+geom_vline(xintercept = 6.5,color="light gray",size = 0.15)+geom_vline(xintercept = 7.5,color="light gray",size = 0.15)+geom_vline(xintercept = 8.5,color="light gray",size = 0.15)+
  geom_bar(stat = "identity", aes(fill=Group), show.legend = F,width = 0.75,position = "dodge",size=0.25) +
  theme_bw() +theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25))+
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),position = "right", 
                     labels = NULL,minor_breaks = NULL)+ theme(panel.grid.major.x = element_blank())+
  scale_x_discrete(labels= NULL)+labs(x=NULL,y=NULL, title =NULL)
Plot_Top2
ggsave(filename = "ALL_T10.EnhAtlas_EnhE.GenicF.enrich.100kb.NO.LABS.pdf",
       path = out,plot = Plot_Top2,device = "pdf",width = 2,height = 1,units = "in")


######
## Bottom regions:
All_df<-as.data.frame(matrix(nrow = 0,ncol = 6)  )

for (i in 1:3){
  Enh_df<-load_obj(f = paste0(path,Enh_Bottom_Names[i]))
  SupEnh_df<-load_obj(f = paste0(path,Sup_Enh_Bottom_Names[i]))
  Genic_df<-load_obj(f = paste0(path,Genic_Names_Bottom[i]))
  
  colnames(Enh_df)<-Enh
  colnames(SupEnh_df)<-SupEnh
  colnames(Genic_df)<-genicF
  Genic_df<-Genic_df[,c(2,1,3:7)]
  
  G_E_df<-cbind(SupEnh_df,Enh_df,Genic_df)
  G_E_df<-Rep_step(G_E_df)
  G_E_df$feature<-factor(x = G_E_df$feature,levels = All_features)
  G_E_df$Group<-Sample_Names[i]
  
  All_df<-rbind(All_df, G_E_df)
  colnames(All_df)<-colnames(G_E_df)
  
}
All_df$Group<-factor(x = All_df$Group,levels = Sample_Names)
All_df$Diff[All_df$Diff >= 6]=6
All_df$Diff[All_df$Diff <= -6]=-6


# BOTTOM regions
Plot_Bottom<-ggplot(data = All_df, aes(x=feature, y=Diff), color=Group) + 
  geom_bar(stat = "identity", aes(fill=Group), show.legend = F,width = 0.75,position = "dodge") +theme_bw() +
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),
                     labels = expression(10^6,10^3, 0,10^3,10^6),
                     minor_breaks = NULL) +
  scale_x_discrete(labels= All_features)+
  theme(plot.title = element_text(size = 10, face = "bold" ,hjust = 0.5),
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 45,hjust = 1),axis.ticks.y.left = element_line(size = 0.25))+
  labs(x=xLab,y=yLab, title = NULL, parse=T)
Plot_Bottom

ggsave(filename = "ALL_B10.EnhAtlas_EnhE.GenicF.enrich.100kb.pdf",
       path = out,plot = Plot_Bottom,device = "pdf",width = 2,height = 1,units = "in")

# BOTTOM regions - plot rotated 90 degrees, no labels:
Plot_Bottom2<-ggplot(data = All_df, aes(x=feature, y=Diff), color=Group) + 
  geom_vline(xintercept = 1.5,color="light gray",size = 0.15)+geom_vline(xintercept = 2.5,color="light gray",size = 0.15)+geom_vline(xintercept = 3.5,color="light gray",size = 0.15)+geom_vline(xintercept = 4.5,color="light gray",size = 0.15)+geom_vline(xintercept = 5.5,color="light gray",size = 0.15)+geom_vline(xintercept = 6.5,color="light gray",size = 0.15)+geom_vline(xintercept = 7.5,color="light gray",size = 0.15)+geom_vline(xintercept = 8.5,color="light gray",size = 0.15)+
  geom_bar(stat = "identity", aes(fill=Group), show.legend = F,width = 0.75,position = "dodge",size=0.25) +
  theme_bw() +theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25))+
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),position = "right", 
                     labels = NULL,minor_breaks = NULL)+ theme(panel.grid.major.x = element_blank())+
  scale_x_discrete(labels= NULL)+labs(x=NULL,y=NULL, title =NULL)
Plot_Bottom2

ggsave(filename = "ALL_B10.EnhAtlas_EnhE.GenicF.enrich.100kb.NO.LABS.pdf",
       path = out,plot = Plot_Bottom2,device = "pdf",width = 2,height = 1,units = "in")




# END

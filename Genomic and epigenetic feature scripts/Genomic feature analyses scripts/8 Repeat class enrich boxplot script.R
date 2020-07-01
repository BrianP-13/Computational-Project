library(ggplot2)
library(IDPmisc)


## Description: Plot with all repeat class enrichment results from all regions of interest.

out<-"output_directory/"


# Common plot labels
Xlab<-"Repeat Class"
Ylab = bquote(atop(bold("Depletion <--> Enrichment"), "-log10(pval)" ))

Sample_Names_Top<-c("Top_64PP","Top_CPD","Top_RankDiff")
Sample_Names_Bottom<-c("Bottom_64PP","Bottom_CPD","Bottom_RankDiff")
Region_Names<-c("6-4PP","CPD","RankDiff")

# Create ordered character list repeat types for plotting
load(file = "file_path/Repeat Classes with subFamily list.RData")
str(Repeat_Class_list)

Repeat_Class_names<-names(Repeat_Class_list)
Repeat_Class_labels<-c("LINE","SINE","LTR","DNA TE","RNA","Satellite","Simple\nrepeat","Low\ncomplexity")


######
## repClass

# Load repClass enrichment results
load(file = "file_path/RepClass_Enrich_Top10perc.100kb.regions.RData")
rClass.Enrich_list


## Top regions
# Format and extract enrichment results from list:
RepClass_Top_df<-as.data.frame(matrix(ncol = 3,nrow = 0))

for (i in 1:3) {
  df<- rClass.Enrich_list[[i]]
  df<-df[c(3,6),]
  df<-log10(df)
  df[df <= 0]<-0
  
  df<-as.data.frame(t(df))
  df$Diff<- df[,1]-df[,2]
  df$Region<-Sample_Names_Top[i]
  
  df$feature<-Repeat_Class_names
  df<-df[,3:5]
  RepClass_Top_df<-rbind(RepClass_Top_df,df,make.row.names=F)
}

# Set limits on max/min enrichment score
RepClass_Top_df$Diff[RepClass_Top_df$Diff >= 6]= 6
RepClass_Top_df$Diff[RepClass_Top_df$Diff <= -6]= -6

RepClass_Top_df$feature<-factor(RepClass_Top_df$feature,levels = Repeat_Class_names)
RepClass_Top_df$Region<-factor(x = RepClass_Top_df$Region, levels = Sample_Names_Top, labels = Region_Names)


## Plot with labels
Plot_Top<-ggplot(data = RepClass_Top_df, aes(x=feature, y=Diff), color=Region) + 
  geom_bar(stat = "identity", aes(fill=Region), show.legend = F,width = 0.75,position = "dodge") +theme_bw() +
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),
                     labels = expression(10^6,10^3, 0,10^3,10^6),
                     minor_breaks = NULL) +
  scale_x_discrete(labels= Repeat_Class_labels)+
  theme(plot.title = element_text(size = 10, face = "bold" ,hjust = 0.5),
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 45,hjust = 1),axis.ticks.y.left = element_line(size = 0.25))+
  labs(x=Xlab,y=Ylab, title = NULL, parse=T)
Plot_Top
ggsave(filename = "ALL_T10.Repeat.Class.enrich.100kb.pdf",
       path = out,plot = Plot_Top,device = "pdf",width = 2,height = 1,units = "in")

# Plots with no labels- plot rotated 90 degrees
Plot_Top2<-ggplot(data = RepClass_Top_df, aes(x=feature, y=Diff), color=Region) + 
  geom_vline(xintercept = 1.5,color="light gray",size = 0.15)+geom_vline(xintercept = 2.5,color="light gray",size = 0.15)+geom_vline(xintercept = 3.5,color="light gray",size = 0.15)+geom_vline(xintercept = 4.5,color="light gray",size = 0.15)+geom_vline(xintercept = 5.5,color="light gray",size = 0.15)+geom_vline(xintercept = 6.5,color="light gray",size = 0.15)+geom_vline(xintercept = 7.5,color="light gray",size = 0.15)+geom_vline(xintercept = 8.5,color="light gray",size = 0.15)+
  geom_bar(stat = "identity", aes(fill=Region), show.legend = F,width = 0.75,position = "dodge",size=0.25) +
  theme_bw() +theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25))+
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),position = "right", 
                     labels = NULL,minor_breaks = NULL)+ theme(panel.grid.major.x = element_blank())+
  scale_x_discrete(labels= NULL)+labs(x=NULL,y=NULL, title =NULL)
Plot_Top2
ggsave(filename = "ALL_T10.Repeat.Class.enrich.100kb.NO.LABS.pdf",
       path = out,plot = Plot_Top2,device = "pdf",width = 2,height = 1,units = "in")


## Bottom regions
# Format and extract enrichment results from list:
RepClass_Bottom_df<-as.data.frame(matrix(ncol = 3,nrow = 0))

for (i in 1:3) {
  df<- rClass.Enrich_list[[i+3]]
  df<-df[c(3,6),]
  df<-log10(df)
  df[df <= 0]<-0
  
  df<-as.data.frame(t(df))
  df$Diff<- df[,1]-df[,2]
  df$Region<-Sample_Names_Bottom[i]
  
  df$feature<-Repeat_Class_names
  df<-df[,3:5]
  RepClass_Bottom_df<-rbind(RepClass_Bottom_df,df,make.row.names=F)
}

# Set limits on max/min enrichment score
RepClass_Bottom_df$Diff[RepClass_Bottom_df$Diff >= 6]= 6
RepClass_Bottom_df$Diff[RepClass_Bottom_df$Diff <= -6]= -6

RepClass_Bottom_df$feature<-factor(RepClass_Bottom_df$feature,levels = Repeat_Class_names)
RepClass_Bottom_df$Region<-factor(x = RepClass_Bottom_df$Region, levels = Sample_Names_Bottom, labels = Region_Names)

## Plot with labels
Plot_Bottom<-ggplot(data = RepClass_Bottom_df, aes(x=feature, y=Diff), color=Region) + 
  geom_bar(stat = "identity", aes(fill=Region), show.legend = F,width = 0.75,position = "dodge") +theme_bw() +
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),
                     labels = expression(10^6,10^3, 0,10^3,10^6),
                     minor_breaks = NULL) +
  scale_x_discrete(labels= Repeat_Class_labels)+
  theme(plot.title = element_text(size = 10, face = "bold" ,hjust = 0.5),
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8, angle = 45,hjust = 1),axis.ticks.y.left = element_line(size = 0.25))+
  labs(x=Xlab,y=Ylab, title = NULL, parse=T)
Plot_Bottom
ggsave(filename = "ALL_B10.Repeat.Class.enrich.100kb.pdf",
       path = out,plot = Plot_Bottom,device = "pdf",width = 2,height = 1,units = "in")

# Plots with no labels- plot rotated 90 degrees
Plot_Bottom2<-ggplot(data = RepClass_Bottom_df, aes(x=feature, y=Diff), color=Region) + 
  geom_vline(xintercept = 1.5,color="light gray",size = 0.15)+geom_vline(xintercept = 2.5,color="light gray",size = 0.15)+geom_vline(xintercept = 3.5,color="light gray",size = 0.15)+geom_vline(xintercept = 4.5,color="light gray",size = 0.15)+geom_vline(xintercept = 5.5,color="light gray",size = 0.15)+geom_vline(xintercept = 6.5,color="light gray",size = 0.15)+geom_vline(xintercept = 7.5,color="light gray",size = 0.15)+geom_vline(xintercept = 8.5,color="light gray",size = 0.15)+
  geom_bar(stat = "identity", aes(fill=Region), show.legend = F,width = 0.75,position = "dodge",size=0.25) +
  theme_bw() +theme(line = element_line(size = 0.15), 
                    panel.border = element_rect(fill = NA,size = 0.25))+
  coord_cartesian(ylim = c(-6,6))+
  scale_y_continuous(breaks = c(-6,-3,0,3,6),position = "right", 
                     labels = NULL,minor_breaks = NULL)+ theme(panel.grid.major.x = element_blank())+
  scale_x_discrete(labels= NULL)+labs(x=NULL,y=NULL, title =NULL)
Plot_Bottom2
ggsave(filename = "ALL_B10.Repeat.Class.enrich.100kb.NO.LABS.pdf",
       path = out,plot = Plot_Bottom2,device = "pdf",width = 2,height = 1,units = "in")



## END ##
library(ggplot2)
library(ggrepel)


# Description: Script measures the mean UV lesion signal per chromosome vs the radial position of each chromosome

out<-"output_directory/"

path<-'file_path/'


# To load completed datatable
RadialPos<-read.csv("file_path/IMR91_chrom_radial_positions.csv")
RadialPos<-RadialPos[,c(2:4)]
colnames(RadialPos)<-c("chrom", "Mean, normalized CT-CN distances", "SD")
chromList<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY")
RadialPos[,1]<-chromList
RadialPos$IP64_mean<-0
RadialPos$CPD_mean<-0
RadialPos$RankDiff_mean<-0
RadialPos<-RadialPos[c(1:22),] # omit row 23, no X chr in this analysis

# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist.RData")
Ranked_diff<-Ranked_diff[which(Ranked_diff$seqnames != "chrX"),]

List<-as.character(unique(Ranked_diff$seqnames))

for (i in 1:22) { 
  df<-subset(x = Ranked_diff,seqnames == List[i])
  
  Mean_64<-mean(x = df$IP64)
  Mean_CPD<-mean(x = df$CPD)
  Mean_Diff<-mean(x = df$Difference)
  
  RadialPos[i,4]<-Mean_64
  RadialPos[i,5]<-Mean_CPD
  RadialPos[i,6]<-Mean_Diff
  
}
RadialPos

#####
## 6-4PP

CorTest<-cor.test(RadialPos[,2], RadialPos[,4])

cor<-CorTest$estimate
corr_coeff<-paste0("r = ",signif(cor, digits = 3))
pval<-CorTest$p.value
pval_coef<- paste0("p = ", print(signif(pval, digits=3)))

xlab<-"Mean 6-4PP signal per chromosome" 
ylab<-"Radial Position"

# 64PP plot with labels
IP64_plot<-ggplot(data = RadialPos, aes(x=RadialPos[,4], y=RadialPos[,2], label=RadialPos[,1]))+ 
  geom_point(size=0.25, shape=19,color="black")+
  labs(x=xlab, y=ylab) + 
  geom_smooth(colour= "red", method=lm)+
  geom_text_repel(aes(label=RadialPos[,1])) +
  annotate(geom="text",label=corr_coeff, size=5,fontface="bold.italic", hjust=0, vjust=-1,
           x = 0.825,y = 0.75) +
  annotate(geom="text",label=pval_coef, size=5, fontface="bold.italic",hjust=0, vjust=0
           ,x = 0.825,y = 0.75) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8))

IP64_plot
ggsave(filename = "64PP_vs_radial_pos.pdf",
       path = out,plot = IP64_plot,device = "pdf",width = 4,height = 4,units = "in")


#####
## CPD
CorTest<-cor.test(RadialPos[,2], RadialPos[,5])

cor<-CorTest$estimate
corr_coeff<-paste0("r = ",signif(cor, digits = 3))
pval<-CorTest$p.value
pval_coef<- paste0("p = ", print(signif(pval, digits=3)))

xlab<-"Mean CPD signal per chromosome" 
ylab<-"Radial Position"

# CPD plot with labels
CPD_plot<-ggplot(data = RadialPos, aes(x=RadialPos[,5], y=RadialPos[,2], label=RadialPos[,1]))+ 
  geom_point(size=0.25, shape=19)+
  labs(x=xlab, y=ylab) + 
  geom_smooth(colour= "red", method=lm)+
  geom_text_repel(aes(label=RadialPos[,1])) +
  annotate(geom="text",label=corr_coeff, size=5,fontface="bold.italic", hjust=0, vjust=-1,
           x = 0.7,y = 0.75) +
  annotate(geom="text",label=pval_coef, size=5, fontface="bold.italic",hjust=0, vjust=0
           ,x = 0.7,y = 0.75) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8))

CPD_plot
ggsave(filename = "CPD_vs_radial_pos.pdf",
       path = out,plot = CPD_plot,device = "pdf",width = 4,height = 4,units = "in")


#####
## RankDiff
CorTest<-cor.test(RadialPos[,2], RadialPos[,6])

cor<-CorTest$estimate
corr_coeff<-paste0("r = ",signif(cor, digits = 3))
pval<-CorTest$p.value
pval_coef<- paste0("p = ", print(signif(pval, digits=3)))

xlab<-"Mean RankDiff signal per chromosome" 

# RankDiff plot with labels
RankDiff_plot<-ggplot(data = RadialPos, aes(x=RadialPos[,6], y=RadialPos[,2], label=RadialPos[,1]))+ 
  geom_point(size=0.25, shape=19)+
  labs(x=xlab, y=ylab) + 
  geom_smooth(colour= "red", method=lm)+
  geom_text_repel(aes(label=RadialPos[,1])) +
  annotate(geom="text",label=corr_coeff, size=5,fontface="bold.italic", hjust=0, vjust=-1,
           x = 0.4,y = 0.75) +
  annotate(geom="text",label=pval_coef, size=5, fontface="bold.italic",hjust=0, vjust=0
           ,x = 0.4,y = 0.75) +
  theme_bw()+theme(line = element_line(size = 0.15),
                   panel.border = element_rect(fill = NA,size = 0.25),
                   axis.title.x = element_text(size = 10),
                   axis.text.x = element_text(size = 8),
                   axis.title.y = element_text(size = 10),
                   axis.text.y = element_text(size = 8))

RankDiff_plot
ggsave(filename = "RankDiff_vs_radial_pos.pdf",
       path = out,plot = RankDiff_plot,device = "pdf",width = 4,height = 4,units = "in")


## END ##

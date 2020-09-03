library(reshape2)
library(ggplot2)


## Description: This script will measure the distance of each TAD from the center of the nucleus (radius) and the rsults will be in boxplot form.
# Will use a custom function to extract the Radial position from the .cmm file
# TAD coordinates are from 3D genomes with TADs with 10% or greater overlap with top 10% regions of interest (64PP,CPD,RankDif) using 100kb bin size

## Strategy:
# Step 1: Will need to parse through the .cmm file containing the bead info (TADs), "rgb" (color) info, and keep only the rows that contain genomic info (no linker regions)
# Step 2: Will need to figure out how to split the character string to separate the relevant coordinate info and add it to a new data frame 

out<-"output_directory/"


# Function for getting TAD coordinate info from CMM files
chrom_3D<-function(cmm_file) { 
  Pattern<-"<marker*"
  param<-c("marker_id","x","y","z","radius","r","g","b","chrID","beadID", "chr","start", "end", "width")
  Chrom3D_TAD<-as.data.frame(matrix(ncol = 14, nrow = 0))
  splitter<-c(" ")
  
  for (i in 1:nrow(cmm_file)) { 
    Row_num<-cmm_file[i,]
    
    if (grepl(pattern =  Pattern, x = Row_num)) {
      splt<-unlist(strsplit(x = as.character(Row_num), split = splitter))
      
      marker<-splt[2]
      marker<-unlist(strsplit(x = marker, split = "="))[2]
      
      x<-splt[3]
      x<-as.numeric(unlist(strsplit(x = x,split = "="))[2])
      
      y<-splt[4]
      y<-as.numeric(unlist(strsplit(x = y,split = "="))[2])
      
      z<-splt[5]
      z<-as.numeric(unlist(strsplit(x = z,split = "="))[2])
      
      radius<-splt[6]
      radius<-unlist(strsplit(x = radius,split = "="))[2]
      
      r<-splt[7]
      r<-as.numeric(unlist(strsplit(x = r,split = "="))[2])
      
      g<-splt[8]
      g<-as.numeric(unlist(strsplit(x = g,split = "="))[2])
      
      b<-splt[9]
      b<-as.numeric(unlist(strsplit(x = b,split = "="))[2])
      
      chrID<-splt[10]
      chrID<-unlist(strsplit(x = chrID,split = "="))[2]
      
      beadID<-splt[11]
      beadID<-unlist(strsplit(x = beadID,split = "="))[2]
      beadID<-unlist(strsplit(x = beadID,split = "/"))[1]
      
      chr<-unlist(strsplit(x = beadID, split = ":"))[1]
      chr<-unlist(strsplit(x = chr, split = "_"))[1]
      
      start<-unlist(strsplit(x = beadID, split = ":"))[2]
      start<-as.numeric(unlist(strsplit(x = start, split = "-"))[1])
      start<-start+1
      
      end<-unlist(strsplit(x = beadID, split = ":"))[2]
      end<-as.numeric(unlist(strsplit(x = end, split = "-"))[2])
      
      width<-end-start+1
      
      Row_param<-c(marker,x,y,z,radius,r,g,b,chrID,beadID,chr,start,end,width)
      Chrom3D_TAD<-rbind.data.frame(Chrom3D_TAD, as.vector(Row_param), stringsAsFactors=F, make.row.names = F)
      #Chrom3D_TAD<-rbind(Chrom3D_TAD, as.vector(Row_param), stringsAsFactors=F) 
      
    } else { 
      next
    }
    
  }
  colnames(Chrom3D_TAD)<-param
  Chrom3D_TAD$x<-as.numeric(Chrom3D_TAD$x)
  Chrom3D_TAD$y<-as.numeric(Chrom3D_TAD$y)
  Chrom3D_TAD$z<-as.numeric(Chrom3D_TAD$z)
  Chrom3D_TAD$r<-as.numeric(Chrom3D_TAD$r)
  Chrom3D_TAD$g<-as.numeric(Chrom3D_TAD$g)
  Chrom3D_TAD$b<-as.numeric(Chrom3D_TAD$b)
  Chrom3D_TAD$start<-as.numeric(Chrom3D_TAD$start)
  Chrom3D_TAD$end<-as.numeric(Chrom3D_TAD$end)
  Chrom3D_TAD$width<-as.numeric(Chrom3D_TAD$width)
  return(Chrom3D_TAD)
  
}


# Function for getting the distance measurements for TADs of interest
TAD_radius<-function(cmm_TAD_info) {
  Top<-cmm_TAD_info[which(cmm_TAD_info$r==1),]
  Top_r<-Top[,c(2:4)]
  Top_r$r<-0
  Bottom<-cmm_TAD_info[which(cmm_TAD_info$b==1),]
  Bottom_r<-Bottom[,c(2:4)]
  Bottom_r$r<-0
  
  for (i in 1:nrow(Top_r)) {
    X<-(Top_r[i,1])^2
    Y<-(Top_r[i,2])^2
    Z<-(Top_r[i,3])^2
    R<-(X+Y+Z)^0.5
    Top_r[i,4]<-R
  }
  
  for (i in 1:nrow(Bottom_r)) {
    X<-(Bottom_r[i,1])^2
    Y<-(Bottom_r[i,2])^2
    Z<-(Bottom_r[i,3])^2
    R<-(X+Y+Z)^0.5
    Bottom_r[i,4]<-R
  }
  
  TOP<-as.data.frame(matrix(nrow = nrow(Top_r), ncol = 2))
  TOP[,1]<-Top_r$r
  TOP[,2]<-"Top_10"
  
  BOTTOM<-as.data.frame(matrix(nrow = nrow(Bottom_r), ncol = 2))
  BOTTOM[,1]<-Bottom_r$r
  BOTTOM[,2]<-"Bottom_10"
  
  TOP_BOTTOM<-rbind(TOP, BOTTOM)
  colnames(TOP_BOTTOM)<-c("Distance","TAD_type")
  TOP_BOTTOM[,2]<-factor(x = TOP_BOTTOM[,2], levels = c('Top_10', 'Bottom_10'))
  #TOP_BOTTOM<-matrix(ncol = 2)
  #TOP_BOTTOM[,1]<-Top_r$r
  #TOP_BOTTOM[,2]<-Bottom_r$r
  #TOP_BOTTOM<-cbind(Top_r$r, Bottom_r$r)
  #TOP_BOTTOM<-data.frame(TOP_BOTTOM)
  #colnames(TOP_BOTTOM)<-c("Top_2SD", "Bottom_2SD")
  return(TOP_BOTTOM)
}


######
## 3D genomes where highlighted TADs (> 10% overlap) represent top/bottom 10% of UV lesions binned to TADs ##
## This also includes 3D genomes generated using signals from top/bottom ranked differences method ##

# Load IMR90 genome with TOP and BOTTOM 10% highlighted TADs
IP64_10<-read.delim(file = "file_path/T_B_10_64PP.100kb_10per_TAD_ovrlp_MERGED.cmm", header = T)
CPD_10<-read.delim(file = "file_path/T_B_10_CPD.100kb_10per_TAD_ovrlp_MERGED.cmm", header = T)
Diff_10<-read.delim(file = "file_path/T_B_10_RankDiff.100kb_10per_TAD_ovrlp_MERGED.cmm", header = T)

# Extract TAD coordinates from .cmm files from highlighted TADs only
IP64_10_TAD<-chrom_3D(cmm_file = IP64_10)
CPD_10_TAD<-chrom_3D(cmm_file = CPD_10)
Diff_10_TAD<-chrom_3D(cmm_file = Diff_10)

# TAD distances from Top and bottom regions
IP64_radius<-TAD_radius(cmm_TAD_info = IP64_10_TAD)
CPD_radius<-TAD_radius(cmm_TAD_info = CPD_10_TAD)
Diff_radius<-TAD_radius(cmm_TAD_info = Diff_10_TAD)


## Wilcox test for significance
Stat_Test<-as.data.frame( matrix(nrow =3 ,ncol = 1) )
colnames(Stat_Test)<-c('p_val_100kb')
rownames(Stat_Test)<-c('IP64','CPD','RankDiff')

IP64_test<-wilcox.test(x = IP64_radius[which(IP64_radius$TAD_type=="Top_10"),1],y = IP64_radius[which(IP64_radius$TAD_type=="Bottom_10"),1], alternative = "two.sided", paired = F, conf.int = T)
Stat_Test[1,]<-IP64_test$p.value

CPD_test<-wilcox.test(x = CPD_radius[which(CPD_radius$TAD_type=="Top_10"),1], y = CPD_radius[which(CPD_radius$TAD_type=="Bottom_10"),1], alternative = "two.sided", paired = F, conf.int = T)
Stat_Test[2,]<-CPD_test$p.value

Diff_test<-wilcox.test(x = Diff_radius[which(Diff_radius$TAD_type=="Top_10"),1],y = Diff_radius[which(Diff_radius$TAD_type=="Bottom_10"),1], alternative = "two.sided", paired = F, conf.int = T)
Stat_Test[3,]<-Diff_test$p.value

write.csv(Stat_Test, file = paste0(out,'T10.10perTAD.ovrlp.wilcoxon_RST.100kb.csv'))


## Plotting script
Ylab<-"Distance from \n center of nucleus"

## 6-4PP plot with labels
IP64_radius$TAD_type<-factor(IP64_radius$TAD_type,levels = c("Bottom_10","Top_10"))

xlab_64<-"6-4PP lesions \n within IMR-90 TADs"
Plot_64<-ggplot(data = IP64_radius, mapping = aes(x = TAD_type, y = Distance, fill=TAD_type)) + 
  geom_boxplot(outlier.shape = NA,size=0.25) + guides(fill=FALSE) +
  xlab(label = xlab_64) + ylab(label = Ylab)+
  annotate(geom = "text", label= IP64_B10_n, x = 1,y=6.5,vjust=-1) +
  annotate(geom = "text", label= IP64_T10_n, x = 2,y=6.5,vjust=-1) +
  
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values = c('sky blue', 'red'))+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
Plot_64
ggsave(filename = "T10_64PP.100kb.TAD.dist.pdf",
       path = out,plot = Plot_64,device = "pdf",width = 1.5,height = 2,units = "in")

#
## CPD plot
CPD_radius$TAD_type<-factor(CPD_radius$TAD_type,levels = c("Bottom_10","Top_10"))
xlab_cpd<-"CPD lesions \n within IMR-90 TADs"

Plot_cpd<-ggplot(data = CPD_radius, mapping = aes(x = TAD_type, y = Distance, fill=TAD_type)) + 
  geom_boxplot(outlier.shape = NA,size=0.25) + guides(fill=FALSE) +
  xlab(label = xlab_cpd) + ylab(label = Ylab)+
  annotate(geom = "text", label= CPD_B10_n, x = 1,y=6.5,vjust=-1) +
  annotate(geom = "text", label= CPD_T10_n, x = 2,y=6.5,vjust=-1) +
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values = c('sky blue', 'red'))+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
Plot_cpd
ggsave(filename = "T10_CPD.100kb.TAD.dist.pdf",
       path = out,plot = Plot_cpd,device = "pdf",width = 1.5,height = 2,units = "in")

#
## RankDiff plots
Diff_radius$TAD_type<-factor(Diff_radius$TAD_type,levels = c("Bottom_10","Top_10"))
xlab_diff<-"TADs with highest ranked \n differences in UV lesion"

Plot_diff<-ggplot(data = Diff_radius, mapping = aes(x = TAD_type, y = Distance, fill=TAD_type)) + 
  geom_boxplot(outlier.shape = NA,size=0.25) + guides(fill=FALSE) +
  xlab(label = xlab_diff)+ylab(label = Ylab)+
  annotate(geom = "text", label= Diff_T10_n, x = 1,y=6.5,vjust=-1) +
  annotate(geom = "text", label= Diff_B10_n, x = 2,y=6.5,vjust=-1) +
  scale_y_continuous(limits = c(0,8))+
  scale_fill_manual(values = c('sky blue', 'red'))+
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     line = element_line(size = 0.25),
                     panel.border = element_rect(fill = NA,size = 0.25),
                     axis.title.x = element_text(size = 10),
                     axis.text.x = element_text(size = 8),
                     axis.title.y = element_text(size = 10),
                     axis.text.y = element_text(size = 8))
Plot_diff
ggsave(filename = "T10_RankDiff.100kb.TAD.dist.pdf",
       path = out,plot = Plot_diff,device = "pdf",width = 1.5,height = 2,units = "in")


## END ##

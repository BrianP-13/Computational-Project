library(GenomicRanges)
library(ggplot2)
library(IDPmisc)
library(reshape2)


out<-"output_directory/"

# Overlap function:
Gene_Ovrlp_fun<-function( Query, Subject,Overlap,Region_name){
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

## Import rank normalized UV dmg dataset and calculate top and bottom 10% regions
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff,keep.extra.columns = T)

# Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))

# Load enhancer dataset 
load(file = 'file_path/Enh_Total_XR_and_mutR_dataset.RData')

# Subset enhancer regions to include only 64PP repair
Enh_XR64<-Enh_XR_mutR[,c(1:7,9)]
Enh_XR64$Log10_mutR<-log10(Enh_XR64$Total_C_to_T_rate)
Enh_XR64$Region<-"All enhancers"

# Subset enhancer regions to include only CPD repair
Enh_XRcpd<-Enh_XR_mutR[,c(1:6,8,9)]
Enh_XRcpd$Log10_mutR<-log10(Enh_XRcpd$Total_C_to_T_rate)
Enh_XRcpd$Region<-"All enhancers"


# Enhancer Mutation stats - All enhancers, with chrX ------------------------

### Stats Table
Enh_XR_mutR_table<-as.data.frame( matrix(nrow = 9,ncol = 5))
rownames(Enh_XR_mutR_table)<-c("All enhancers","Enhancers without 6-4PP repair","Enhancers with 6-4PP repair","Percent without 6-4PP repair","Percent with 6-4PP repair",
                               "Enhancers without CPD repair","Enhancers with CPD repair","Percent without CPD repair","Percent with CPD repair")
colnames(Enh_XR_mutR_table)<-c("All enhancers","Enhancers with mutation rate > 0","Enhancers with no mutation rate","Percent mutated", "Percent not mutated")

# Proportion of enhancers with 64 and CPD repair
Enh_XR_mutR_table[1,1]<-nrow(Enh_XR_mutR)
Enh_XR_mutR_table[2,1]<-nrow(Enh_XR_mutR[which(Enh_XR_mutR$Total_XR64_strand_mean_repair==0),]) # Number of enhancers without 64PP repair
Enh_XR_mutR_table[3,1]<- nrow(Enh_XR_mutR[which(Enh_XR_mutR$Total_XR64_strand_mean_repair!=0),]) # Number of enhancers with 64PP repair
Enh_XR_mutR_table[4,1]<-(Enh_XR_mutR_table[2,1]/nrow(Enh_XR_mutR))*100 # % of enhancers without 64PP repair
Enh_XR_mutR_table[5,1]<-(Enh_XR_mutR_table[3,1]/nrow(Enh_XR_mutR))*100 # % of enhancers with 64PP repair
Enh_XR_mutR_table[6,1]<-nrow(Enh_XR_mutR[which(Enh_XR_mutR$Total_XRCPD_strand_mean_repair==0),]) # Number of enhancers without CPD repair
Enh_XR_mutR_table[7,1]<- nrow(Enh_XR_mutR[which(Enh_XR_mutR$Total_XRCPD_strand_mean_repair!=0),]) # Number of enhancers with CPD repair
Enh_XR_mutR_table[8,1]<-(Enh_XR_mutR_table[6,1]/nrow(Enh_XR_mutR))*100 # % of enhancers without CPD repair
Enh_XR_mutR_table[9,1]<-(Enh_XR_mutR_table[7,1]/nrow(Enh_XR_mutR))*100 # % of enhancers with CPD repair

# subset ehnancers with C>T mutations
Enh_XR_w_MutR<-subset(Enh_XR_mutR,Total_C_to_T_rate > 0)

Enh_XR_mutR_table[1,2]<-nrow(Enh_XR_w_MutR) # All Enhancers with C>T mutations
Enh_XR_mutR_table[2,2]<-nrow(Enh_XR_w_MutR[which(Enh_XR_w_MutR$Total_XR64_strand_mean_repair==0),]) # Mutated enhancers without 64PP repair
Enh_XR_mutR_table[3,2]<-nrow(Enh_XR_w_MutR[which(Enh_XR_w_MutR$Total_XR64_strand_mean_repair!=0),]) # Mutated enhancers with 64PP repair
Enh_XR_mutR_table[4,2]<-(Enh_XR_mutR_table[2,2]/nrow(Enh_XR_w_MutR))*100 # Percent with no 64PP repair
Enh_XR_mutR_table[5,2]<-(Enh_XR_mutR_table[3,2]/nrow(Enh_XR_w_MutR))*100 # Percent with 64PP repair
Enh_XR_mutR_table[6,2]<-nrow(Enh_XR_w_MutR[which(Enh_XR_w_MutR$Total_XRCPD_strand_mean_repair==0),]) # Mutated enhancers without CPD repair
Enh_XR_mutR_table[7,2]<-nrow(Enh_XR_w_MutR[which(Enh_XR_w_MutR$Total_XRCPD_strand_mean_repair!=0),]) # Mutated enhancers with CPD repair
Enh_XR_mutR_table[8,2]<-(Enh_XR_mutR_table[6,2]/nrow(Enh_XR_w_MutR))*100 # Percent with no CPD repair
Enh_XR_mutR_table[9,2]<-(Enh_XR_mutR_table[7,2]/nrow(Enh_XR_w_MutR))*100 # Percent with CPD repair

# subset ehnancers without C>T mutations
Enh_XR_No_MutR<-subset(Enh_XR_mutR,Total_C_to_T_rate == 0)

Enh_XR_mutR_table[1,3]<-nrow(Enh_XR_No_MutR) # All Enhancers with C>T mutations
Enh_XR_mutR_table[2,3]<-nrow(Enh_XR_No_MutR[which(Enh_XR_No_MutR$Total_XR64_strand_mean_repair==0),]) # Mutated enhancers without 64PP repair
Enh_XR_mutR_table[3,3]<-nrow(Enh_XR_No_MutR[which(Enh_XR_No_MutR$Total_XR64_strand_mean_repair!=0),]) # Mutated enhancers with 64PP repair
Enh_XR_mutR_table[4,3]<-(Enh_XR_mutR_table[2,3]/nrow(Enh_XR_No_MutR))*100 # Percent with no 64PP repair
Enh_XR_mutR_table[5,3]<-(Enh_XR_mutR_table[3,3]/nrow(Enh_XR_No_MutR))*100 # Percent with 64PP repair
Enh_XR_mutR_table[6,3]<-nrow(Enh_XR_No_MutR[which(Enh_XR_No_MutR$Total_XRCPD_strand_mean_repair==0),]) # Mutated enhancers without CPD repair
Enh_XR_mutR_table[7,3]<-nrow(Enh_XR_No_MutR[which(Enh_XR_No_MutR$Total_XRCPD_strand_mean_repair!=0),]) # Mutated enhancers with CPD repair
Enh_XR_mutR_table[8,3]<-(Enh_XR_mutR_table[6,3]/nrow(Enh_XR_No_MutR))*100 # Percent with no CPD repair
Enh_XR_mutR_table[9,3]<-(Enh_XR_mutR_table[7,3]/nrow(Enh_XR_No_MutR))*100 # Percent with CPD repair

## Remaining stats: Percentage of enhancers (all, without and with repair) that have C>T mutations (or no mutations)
Enh_XR_mutR_table[1,4]<-(Enh_XR_mutR_table[1,2]/Enh_XR_mutR_table[1,1])*100 # Percent of all enhancers mutated
Enh_XR_mutR_table[1,5]<-(Enh_XR_mutR_table[1,3]/Enh_XR_mutR_table[1,1])*100 # Percent of all enhancers not mutated

Enh_XR_mutR_table[2,4]<-(Enh_XR_mutR_table[2,2]/Enh_XR_mutR_table[2,1])*100 # Percent of enhancers with NO 64PP repair mutated
Enh_XR_mutR_table[2,5]<-(Enh_XR_mutR_table[2,3]/Enh_XR_mutR_table[2,1])*100 # Percent of enhancers with NO 64PP repair not mutated

Enh_XR_mutR_table[3,4]<-(Enh_XR_mutR_table[3,2]/Enh_XR_mutR_table[3,1])*100 # Percent of enhancers with 64PP repair mutated
Enh_XR_mutR_table[3,5]<-(Enh_XR_mutR_table[3,3]/Enh_XR_mutR_table[3,1])*100 # Percent of enhancers with 64PP repair not mutated

Enh_XR_mutR_table[6,4]<-(Enh_XR_mutR_table[6,2]/Enh_XR_mutR_table[6,1])*100 # Percent of enhancers with NO CPD repair mutated
Enh_XR_mutR_table[6,5]<-(Enh_XR_mutR_table[6,3]/Enh_XR_mutR_table[6,1])*100 # Percent of enhancers with NO CPD repair not mutated

Enh_XR_mutR_table[7,4]<-(Enh_XR_mutR_table[7,2]/Enh_XR_mutR_table[7,1])*100 # Percent of enhancers with CPD repair mutated
Enh_XR_mutR_table[7,5]<-(Enh_XR_mutR_table[7,3]/Enh_XR_mutR_table[7,1])*100 # Percent of enhancers with CPD repair not mutated

write.csv(Enh_XR_mutR_table, file = paste0(out,"All enhancers - Excision repair and mutation stats - w ChrX.csv"))


# 6-4PP repair analysis - Top 10% 64PP ------------------------------------

## TABLE with enhancer stats in top 64PP regions
T10_Enh_64PP_XR_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(T10_Enh_64PP_XR_mutR_table)<-c("All Enh in Top 64PP","Enhancers without 6-4PP repair","Enhancers with 6-4PP repair","Percent  without repair","Percent  with repair")
colnames(T10_Enh_64PP_XR_mutR_table)<-c("All Enh in Top 64PP","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Top 10%
XR64_T10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_T10,Overlap = 50,Region_name = "Top 10% 64PP")
XR64_T10_Enh$overlap<-NULL
XR64_T10_Enh$percent_overlap<-NULL
XR64_T10_Enh$Region<-NULL
XR64_T10_Enh<-as.data.frame(XR64_T10_Enh)

# Stats
T10_Enh_64PP_XR_mutR_table[1,1]<-nrow(XR64_T10_Enh) # All enhancers in top region
T10_Enh_64PP_XR_mutR_table[1,2]<-nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
T10_Enh_64PP_XR_mutR_table[1,3]<-nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
T10_Enh_64PP_XR_mutR_table[1,4]<-(nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_C_to_T_rate!=0),])/nrow(XR64_T10_Enh))*100 # Percent of all enhancers mutated
T10_Enh_64PP_XR_mutR_table[1,5]<-(nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_C_to_T_rate==0),])/nrow(XR64_T10_Enh))*100 # Percent of all enhancers NOT mutations

T10_Enh_64PP_XR_mutR_table[4,1]<-(nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair==0),])/nrow(XR64_T10_Enh))*100 # Percent enhancers without repair
T10_Enh_64PP_XR_mutR_table[5,1]<-(nrow(XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair!=0),])/nrow(XR64_T10_Enh))*100 # Percent enhancers with repair

# Enhancers in top 10% 64PP regions - no 64PP repair
XR64_T10_Enh_no_repair<-XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair==0),]
XR64_T10_Enh_no_repair$Group<-"Top 10% 64PP - no repair"

T10_Enh_64PP_XR_mutR_table[2,1]<-nrow(XR64_T10_Enh_no_repair) # All enhancers with no repair in top 64 regions
T10_Enh_64PP_XR_mutR_table[2,2]<-nrow(XR64_T10_Enh_no_repair[which(XR64_T10_Enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in top 64 with C>T mutations
T10_Enh_64PP_XR_mutR_table[2,3]<-nrow(XR64_T10_Enh_no_repair[which(XR64_T10_Enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in top 64 with NO C>T mutations
T10_Enh_64PP_XR_mutR_table[2,4]<-(nrow(XR64_T10_Enh_no_repair[which(XR64_T10_Enh_no_repair$Total_C_to_T_rate!=0),])/nrow(XR64_T10_Enh_no_repair))*100 # Percent of enh with no repair in top 64 with C>T mutations
T10_Enh_64PP_XR_mutR_table[2,5]<-(nrow(XR64_T10_Enh_no_repair[which(XR64_T10_Enh_no_repair$Total_C_to_T_rate==0),])/nrow(XR64_T10_Enh_no_repair))*100 # Percent of enh with no repair in top 64 with NO C>T mutations

# Enhancers in top 10% 64PP regions - with 64PP repair
XR64_T10_Enh_w_repair<-XR64_T10_Enh[which(XR64_T10_Enh$Total_XR64_strand_mean_repair!=0),]
XR64_T10_Enh_w_repair$Group<-"Top 10% 64PP - 64PP repair"

T10_Enh_64PP_XR_mutR_table[3,1]<-nrow(XR64_T10_Enh_w_repair) # All enhancers with 64PP repair in top 64 regions
T10_Enh_64PP_XR_mutR_table[3,2]<-nrow(XR64_T10_Enh_w_repair[which(XR64_T10_Enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with 64PP repair in top 64 with C>T mutations
T10_Enh_64PP_XR_mutR_table[3,3]<-nrow(XR64_T10_Enh_w_repair[which(XR64_T10_Enh_w_repair$Total_C_to_T_rate==0),]) # Enh with 64PP repair in top 64 with NO C>T mutations
T10_Enh_64PP_XR_mutR_table[3,4]<-(nrow(XR64_T10_Enh_w_repair[which(XR64_T10_Enh_w_repair$Total_C_to_T_rate!=0),])/nrow(XR64_T10_Enh_w_repair))*100 # Percent of enh with 64PP repair in top 64 with C>T mutations
T10_Enh_64PP_XR_mutR_table[3,5]<-(nrow(XR64_T10_Enh_w_repair[which(XR64_T10_Enh_w_repair$Total_C_to_T_rate==0),])/nrow(XR64_T10_Enh_w_repair))*100 # Percent of enh with 64PP repair in top 64 with NO C>T mutations

# Remaining stats
T10_Enh_64PP_XR_mutR_table[4,2]<-(T10_Enh_64PP_XR_mutR_table[2,2]/T10_Enh_64PP_XR_mutR_table[1,2])*100 # Percent mutated enhancers without 64PP repair
T10_Enh_64PP_XR_mutR_table[5,2]<-(T10_Enh_64PP_XR_mutR_table[3,2]/T10_Enh_64PP_XR_mutR_table[1,2])*100 # Percent mutated enhancers with 64PP repair
T10_Enh_64PP_XR_mutR_table[4,3]<-(T10_Enh_64PP_XR_mutR_table[2,3]/T10_Enh_64PP_XR_mutR_table[1,3])*100 # Percent enhancers not mutated without 64PP repair
T10_Enh_64PP_XR_mutR_table[5,3]<-(T10_Enh_64PP_XR_mutR_table[3,3]/T10_Enh_64PP_XR_mutR_table[1,3])*100 # Percent enhancers not mutated with 64PP repair

write.csv(T10_Enh_64PP_XR_mutR_table, file = paste0(out,"T10 64PP enhancers - XR64 repair and mutation stats - w ChrX.csv"))


# 6-4PP repair analysis - Bottom 10% 6-4PP --------------------------------

## Bottom 10%
XR64_B10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = IP64_B10,Overlap = 50,Region_name = "Bottom 10% 64PP")
XR64_B10_Enh$overlap<-NULL
XR64_B10_Enh$percent_overlap<-NULL
XR64_B10_Enh$Region<-NULL
XR64_B10_Enh<-as.data.frame(XR64_B10_Enh)

# TABLE with enhancer stats in bottom 64PP regions
B10_Enh_64PP_XR_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(B10_Enh_64PP_XR_mutR_table)<-c("All Enh in Bottom 64PP","Enhancers without 6-4PP repair","Enhancers with 6-4PP repair","Percent  without repair","Percent  with repair")
colnames(B10_Enh_64PP_XR_mutR_table)<-c("All Enh in Bottom 64PP","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

# Stats
B10_Enh_64PP_XR_mutR_table[1,1]<-nrow(XR64_B10_Enh) # All enhancers in top region
B10_Enh_64PP_XR_mutR_table[1,2]<-nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
B10_Enh_64PP_XR_mutR_table[1,3]<-nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
B10_Enh_64PP_XR_mutR_table[1,4]<-(nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_C_to_T_rate!=0),])/nrow(XR64_B10_Enh))*100 # Percent of all enhancers mutated
B10_Enh_64PP_XR_mutR_table[1,5]<-(nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_C_to_T_rate==0),])/nrow(XR64_B10_Enh))*100 # Percent of all enhancers NOT mutations

B10_Enh_64PP_XR_mutR_table[4,1]<-(nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair==0),])/nrow(XR64_B10_Enh))*100 # Percent enhancers without repair
B10_Enh_64PP_XR_mutR_table[5,1]<-(nrow(XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair!=0),])/nrow(XR64_B10_Enh))*100 # Percent enhancers with repair

# Enhancers in bottom 10% 64PP regions - no 64PP repair
XR64_B10_Enh_no_repair<-XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair==0),]
XR64_B10_Enh_no_repair$Group<-"Bottom 10% 64PP - no repair"

B10_Enh_64PP_XR_mutR_table[2,1]<-nrow(XR64_B10_Enh_no_repair) # All enhancers with no repair in bottom 64 regions
B10_Enh_64PP_XR_mutR_table[2,2]<-nrow(XR64_B10_Enh_no_repair[which(XR64_B10_Enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in bottom 64 with C>T mutations
B10_Enh_64PP_XR_mutR_table[2,3]<-nrow(XR64_B10_Enh_no_repair[which(XR64_B10_Enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in bottom 64 with NO C>T mutations
B10_Enh_64PP_XR_mutR_table[2,4]<-(nrow(XR64_B10_Enh_no_repair[which(XR64_B10_Enh_no_repair$Total_C_to_T_rate!=0),])/nrow(XR64_B10_Enh_no_repair))*100 # Percent of enh with no repair in bottom 64 with C>T mutations
B10_Enh_64PP_XR_mutR_table[2,5]<-(nrow(XR64_B10_Enh_no_repair[which(XR64_B10_Enh_no_repair$Total_C_to_T_rate==0),])/nrow(XR64_B10_Enh_no_repair))*100 # Percent of enh with no repair in bottom 64 with NO C>T mutations

# Enhancers in bottom 10% 64PP regions - with 64PP repair
XR64_B10_Enh_w_repair<-XR64_B10_Enh[which(XR64_B10_Enh$Total_XR64_strand_mean_repair!=0),]
XR64_B10_Enh_w_repair$Group<-"Bottom 10% 64PP - 64PP repair"

B10_Enh_64PP_XR_mutR_table[3,1]<-nrow(XR64_B10_Enh_w_repair) # All enhancers with 64PP repair in bottom 64 regions
B10_Enh_64PP_XR_mutR_table[3,2]<-nrow(XR64_B10_Enh_w_repair[which(XR64_B10_Enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with 64PP repair in bottom 64 with C>T mutations
B10_Enh_64PP_XR_mutR_table[3,3]<-nrow(XR64_B10_Enh_w_repair[which(XR64_B10_Enh_w_repair$Total_C_to_T_rate==0),]) # Enh with 64PP repair in bottom 64 with NO C>T mutations
B10_Enh_64PP_XR_mutR_table[3,4]<-(nrow(XR64_B10_Enh_w_repair[which(XR64_B10_Enh_w_repair$Total_C_to_T_rate!=0),])/nrow(XR64_B10_Enh_w_repair))*100 # Percent of enh with 64PP repair in bottom 64 with C>T mutations
B10_Enh_64PP_XR_mutR_table[3,5]<-(nrow(XR64_B10_Enh_w_repair[which(XR64_B10_Enh_w_repair$Total_C_to_T_rate==0),])/nrow(XR64_B10_Enh_w_repair))*100 # Percent of enh with 64PP repair in bottom 64 with NO C>T mutations

# Remaining stats
B10_Enh_64PP_XR_mutR_table[4,2]<-(B10_Enh_64PP_XR_mutR_table[2,2]/B10_Enh_64PP_XR_mutR_table[1,2])*100 # Percent mutated enhancers without 64PP repair
B10_Enh_64PP_XR_mutR_table[5,2]<-(B10_Enh_64PP_XR_mutR_table[3,2]/B10_Enh_64PP_XR_mutR_table[1,2])*100 # Percent mutated enhancers with 64PP repair
B10_Enh_64PP_XR_mutR_table[4,3]<-(B10_Enh_64PP_XR_mutR_table[2,3]/B10_Enh_64PP_XR_mutR_table[1,3])*100 # Percent enhancers not mutated without 64PP repair
B10_Enh_64PP_XR_mutR_table[5,3]<-(B10_Enh_64PP_XR_mutR_table[3,3]/B10_Enh_64PP_XR_mutR_table[1,3])*100 # Percent enhancers not mutated with 64PP repair

write.csv(B10_Enh_64PP_XR_mutR_table, file = paste0(out,"B10 64PP enhancers - XR64 repair and mutation stats - w ChrX.csv"))


# CPD repair analysis - Top 10% -------------------------------------------

## TABLE with enhancer stats in top CPD regions
T10_Enh_CPD_XR_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(T10_Enh_CPD_XR_mutR_table)<-c("All Enh in Top CPD","Enhancers without CPD repair","Enhancers with CPD repair","Percent  without repair","Percent  with repair")
colnames(T10_Enh_CPD_XR_mutR_table)<-c("All Enh in Top CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Top 10%
XRcpd_T10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = CPD_T10,Overlap = 50,Region_name = "Top 10% CPD")
XRcpd_T10_Enh$overlap<-NULL
XRcpd_T10_Enh$percent_overlap<-NULL
XRcpd_T10_Enh$Region<-NULL
XRcpd_T10_Enh<-as.data.frame(XRcpd_T10_Enh)

# Stats
T10_Enh_CPD_XR_mutR_table[1,1]<-nrow(XRcpd_T10_Enh) # All enhancers in top region
T10_Enh_CPD_XR_mutR_table[1,2]<-nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
T10_Enh_CPD_XR_mutR_table[1,3]<-nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
T10_Enh_CPD_XR_mutR_table[1,4]<-(nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_C_to_T_rate!=0),])/nrow(XRcpd_T10_Enh))*100 # Percent of all enhancers mutated
T10_Enh_CPD_XR_mutR_table[1,5]<-(nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_C_to_T_rate==0),])/nrow(XRcpd_T10_Enh))*100 # Percent of all enhancers NOT mutations

T10_Enh_CPD_XR_mutR_table[4,1]<-(nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair==0),])/nrow(XRcpd_T10_Enh))*100 # Percent enhancers without repair
T10_Enh_CPD_XR_mutR_table[5,1]<-(nrow(XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair!=0),])/nrow(XRcpd_T10_Enh))*100 # Percent enhancers with repair

# Enhancers in top 10% CPD regions - no CPD repair
XRcpd_T10_Enh_no_repair<-XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair==0),]
XRcpd_T10_Enh_no_repair$Group<-"Top 10% CPD - no repair"

T10_Enh_CPD_XR_mutR_table[2,1]<-nrow(XRcpd_T10_Enh_no_repair) # All enhancers with no repair in top CPD regions
T10_Enh_CPD_XR_mutR_table[2,2]<-nrow(XRcpd_T10_Enh_no_repair[which(XRcpd_T10_Enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in top CPD with C>T mutations
T10_Enh_CPD_XR_mutR_table[2,3]<-nrow(XRcpd_T10_Enh_no_repair[which(XRcpd_T10_Enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in top CPD with NO C>T mutations
T10_Enh_CPD_XR_mutR_table[2,4]<-(nrow(XRcpd_T10_Enh_no_repair[which(XRcpd_T10_Enh_no_repair$Total_C_to_T_rate!=0),])/nrow(XRcpd_T10_Enh_no_repair))*100 # Percent of enh with no repair in top CPD with C>T mutations
T10_Enh_CPD_XR_mutR_table[2,5]<-(nrow(XRcpd_T10_Enh_no_repair[which(XRcpd_T10_Enh_no_repair$Total_C_to_T_rate==0),])/nrow(XRcpd_T10_Enh_no_repair))*100 # Percent of enh with no repair in top CPD with NO C>T mutations

# Enhancers in top 10% CPD regions - with CPD repair
XRcpd_T10_Enh_w_repair<-XRcpd_T10_Enh[which(XRcpd_T10_Enh$Total_XRCPD_strand_mean_repair!=0),]
XRcpd_T10_Enh_w_repair$Group<-"Top 10% CPD - CPD repair"

T10_Enh_CPD_XR_mutR_table[3,1]<-nrow(XRcpd_T10_Enh_w_repair) # All enhancers with CPD repair in top CPD regions
T10_Enh_CPD_XR_mutR_table[3,2]<-nrow(XRcpd_T10_Enh_w_repair[which(XRcpd_T10_Enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with CPD repair in top CPD with C>T mutations
T10_Enh_CPD_XR_mutR_table[3,3]<-nrow(XRcpd_T10_Enh_w_repair[which(XRcpd_T10_Enh_w_repair$Total_C_to_T_rate==0),]) # Enh with CPD repair in top CPD with NO C>T mutations
T10_Enh_CPD_XR_mutR_table[3,4]<-(nrow(XRcpd_T10_Enh_w_repair[which(XRcpd_T10_Enh_w_repair$Total_C_to_T_rate!=0),])/nrow(XRcpd_T10_Enh_w_repair))*100 # Percent of enh with CPD repair in top CPD with C>T mutations
T10_Enh_CPD_XR_mutR_table[3,5]<-(nrow(XRcpd_T10_Enh_w_repair[which(XRcpd_T10_Enh_w_repair$Total_C_to_T_rate==0),])/nrow(XRcpd_T10_Enh_w_repair))*100 # Percent of enh with CPD repair in top CPD with NO C>T mutations

# Remaining stats
T10_Enh_CPD_XR_mutR_table[4,2]<-(T10_Enh_CPD_XR_mutR_table[2,2]/T10_Enh_CPD_XR_mutR_table[1,2])*100 # Percent mutated enhancers without CPD repair
T10_Enh_CPD_XR_mutR_table[5,2]<-(T10_Enh_CPD_XR_mutR_table[3,2]/T10_Enh_CPD_XR_mutR_table[1,2])*100 # Percent mutated enhancers with CPD repair
T10_Enh_CPD_XR_mutR_table[4,3]<-(T10_Enh_CPD_XR_mutR_table[2,3]/T10_Enh_CPD_XR_mutR_table[1,3])*100 # Percent enhancers not mutated without CPD repair
T10_Enh_CPD_XR_mutR_table[5,3]<-(T10_Enh_CPD_XR_mutR_table[3,3]/T10_Enh_CPD_XR_mutR_table[1,3])*100 # Percent enhancers not mutated with CPD repair

write.csv(T10_Enh_CPD_XR_mutR_table, file = paste0(out,"T10 CPD enhancers - XRcpd repair and mutation stats - w ChrX.csv"))


# CPD repair analysis - Bottom 10% ----------------------------------------

## Bottom 10%
XRcpd_B10_Enh<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = CPD_B10,Overlap = 50,Region_name = "Bottom 10% CPD")
XRcpd_B10_Enh$overlap<-NULL
XRcpd_B10_Enh$percent_overlap<-NULL
XRcpd_B10_Enh$Region<-NULL
XRcpd_B10_Enh<-as.data.frame(XRcpd_B10_Enh)

# TABLE with enhancer stats in bottom CPD regions
B10_Enh_CPD_XR_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(B10_Enh_CPD_XR_mutR_table)<-c("All Enh in Bottom CPD","Enhancers without CPD repair","Enhancers with CPD repair","Percent  without repair","Percent  with repair")
colnames(B10_Enh_CPD_XR_mutR_table)<-c("All Enh in Bottom CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

# Stats
B10_Enh_CPD_XR_mutR_table[1,1]<-nrow(XRcpd_B10_Enh) # All enhancers in top region
B10_Enh_CPD_XR_mutR_table[1,2]<-nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
B10_Enh_CPD_XR_mutR_table[1,3]<-nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
B10_Enh_CPD_XR_mutR_table[1,4]<-(nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_C_to_T_rate!=0),])/nrow(XRcpd_B10_Enh))*100 # Percent of all enhancers mutated
B10_Enh_CPD_XR_mutR_table[1,5]<-(nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_C_to_T_rate==0),])/nrow(XRcpd_B10_Enh))*100 # Percent of all enhancers NOT mutations

B10_Enh_CPD_XR_mutR_table[4,1]<-(nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair==0),])/nrow(XRcpd_B10_Enh))*100 # Percent enhancers without repair
B10_Enh_CPD_XR_mutR_table[5,1]<-(nrow(XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair!=0),])/nrow(XRcpd_B10_Enh))*100 # Percent enhancers with repair

# Enhancers in bottom 10% CPD regions - no CPD repair
XRcpd_B10_Enh_no_repair<-XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair==0),]
XRcpd_B10_Enh_no_repair$Group<-"Bottom 10% CPD - no repair"

B10_Enh_CPD_XR_mutR_table[2,1]<-nrow(XRcpd_B10_Enh_no_repair) # All enhancers with no repair in bottom CPD regions
B10_Enh_CPD_XR_mutR_table[2,2]<-nrow(XRcpd_B10_Enh_no_repair[which(XRcpd_B10_Enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in bottom CPD with C>T mutations
B10_Enh_CPD_XR_mutR_table[2,3]<-nrow(XRcpd_B10_Enh_no_repair[which(XRcpd_B10_Enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in bottom CPD with NO C>T mutations
B10_Enh_CPD_XR_mutR_table[2,4]<-(nrow(XRcpd_B10_Enh_no_repair[which(XRcpd_B10_Enh_no_repair$Total_C_to_T_rate!=0),])/nrow(XRcpd_B10_Enh_no_repair))*100 # Percent of enh with no repair in bottom CPD with C>T mutations
B10_Enh_CPD_XR_mutR_table[2,5]<-(nrow(XRcpd_B10_Enh_no_repair[which(XRcpd_B10_Enh_no_repair$Total_C_to_T_rate==0),])/nrow(XRcpd_B10_Enh_no_repair))*100 # Percent of enh with no repair in bottom CPD with NO C>T mutations

# Enhancers in bottom 10% CPD regions - with CPD repair
XRcpd_B10_Enh_w_repair<-XRcpd_B10_Enh[which(XRcpd_B10_Enh$Total_XRCPD_strand_mean_repair!=0),]
XRcpd_B10_Enh_w_repair$Group<-"Bottom 10% CPD - CPD repair"

B10_Enh_CPD_XR_mutR_table[3,1]<-nrow(XRcpd_B10_Enh_w_repair) # All enhancers with CPD repair in bottom CPD regions
B10_Enh_CPD_XR_mutR_table[3,2]<-nrow(XRcpd_B10_Enh_w_repair[which(XRcpd_B10_Enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with CPD repair in bottom CPD with C>T mutations
B10_Enh_CPD_XR_mutR_table[3,3]<-nrow(XRcpd_B10_Enh_w_repair[which(XRcpd_B10_Enh_w_repair$Total_C_to_T_rate==0),]) # Enh with CPD repair in bottom CPD with NO C>T mutations
B10_Enh_CPD_XR_mutR_table[3,4]<-(nrow(XRcpd_B10_Enh_w_repair[which(XRcpd_B10_Enh_w_repair$Total_C_to_T_rate!=0),])/nrow(XRcpd_B10_Enh_w_repair))*100 # Percent of enh with CPD repair in bottom CPD with C>T mutations
B10_Enh_CPD_XR_mutR_table[3,5]<-(nrow(XRcpd_B10_Enh_w_repair[which(XRcpd_B10_Enh_w_repair$Total_C_to_T_rate==0),])/nrow(XRcpd_B10_Enh_w_repair))*100 # Percent of enh with CPD repair in bottom CPD with NO C>T mutations

# Remaining stats
B10_Enh_CPD_XR_mutR_table[4,2]<-(B10_Enh_CPD_XR_mutR_table[2,2]/B10_Enh_CPD_XR_mutR_table[1,2])*100 # Percent mutated enhancers without CPD repair
B10_Enh_CPD_XR_mutR_table[5,2]<-(B10_Enh_CPD_XR_mutR_table[3,2]/B10_Enh_CPD_XR_mutR_table[1,2])*100 # Percent mutated enhancers with CPD repair
B10_Enh_CPD_XR_mutR_table[4,3]<-(B10_Enh_CPD_XR_mutR_table[2,3]/B10_Enh_CPD_XR_mutR_table[1,3])*100 # Percent enhancers not mutated without CPD repair
B10_Enh_CPD_XR_mutR_table[5,3]<-(B10_Enh_CPD_XR_mutR_table[3,3]/B10_Enh_CPD_XR_mutR_table[1,3])*100 # Percent enhancers not mutated with CPD repair

write.csv(B10_Enh_CPD_XR_mutR_table, file = paste0(out,"B10 CPD enhancers - XRcpd repair and mutation stats - w ChrX.csv"))


# 6-4PP repair analysis - Top 10% 6-4PP-CPD --------------------------------

## TABLE with enhancer stats in top 6-4PP-CPD regions
T10_Enh_Diff_XR64_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(T10_Enh_Diff_XR64_mutR_table)<-c("All Enh in Top 6-4PP-CPD","Enhancers without 6-4PP repair","Enhancers with 6-4PP repair","Percent  without repair","Percent  with repair")
colnames(T10_Enh_Diff_XR64_mutR_table)<-c("All Enh in Top 6-4PP-CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Top 10%
T10_Diff_enh_XR64<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% 6-4PP-CPD")
T10_Diff_enh_XR64$overlap<-NULL
T10_Diff_enh_XR64$percent_overlap<-NULL
T10_Diff_enh_XR64$Region<-NULL
T10_Diff_enh_XR64<-as.data.frame(T10_Diff_enh_XR64)

# Stats
T10_Enh_Diff_XR64_mutR_table[1,1]<-nrow(T10_Diff_enh_XR64) # All enhancers in top region
T10_Enh_Diff_XR64_mutR_table[1,2]<-nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
T10_Enh_Diff_XR64_mutR_table[1,3]<-nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
T10_Enh_Diff_XR64_mutR_table[1,4]<-(nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_C_to_T_rate!=0),])/nrow(T10_Diff_enh_XR64))*100 # Percent of all enhancers mutated
T10_Enh_Diff_XR64_mutR_table[1,5]<-(nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_C_to_T_rate==0),])/nrow(T10_Diff_enh_XR64))*100 # Percent of all enhancers NOT mutations

T10_Enh_Diff_XR64_mutR_table[4,1]<-(nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_XR64_strand_mean_repair==0),])/nrow(T10_Diff_enh_XR64))*100 # Percent enhancers without repair
T10_Enh_Diff_XR64_mutR_table[5,1]<-(nrow(T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_XR64_strand_mean_repair!=0),])/nrow(T10_Diff_enh_XR64))*100 # Percent enhancers with repair

# Enhancers in top 10% 6-4PP-CPD regions - no 64PP repair
T10_Diff_enh_no_repair<-T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_XR64_strand_mean_repair==0),]
T10_Diff_enh_no_repair$Group<-"Top 10% 6-4PP-CPD - no repair"

T10_Enh_Diff_XR64_mutR_table[2,1]<-nrow(T10_Diff_enh_no_repair) # All enhancers with no repair in top 6-4PP-CPD regions
T10_Enh_Diff_XR64_mutR_table[2,2]<-nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XR64_mutR_table[2,3]<-nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in top 6-4PP-CPD with NO C>T mutations
T10_Enh_Diff_XR64_mutR_table[2,4]<-(nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate!=0),])/nrow(T10_Diff_enh_no_repair))*100 # Percent of enh with no repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XR64_mutR_table[2,5]<-(nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate==0),])/nrow(T10_Diff_enh_no_repair))*100 # Percent of enh with no repair in top 6-4PP-CPD with NO C>T mutations

# Enhancers in top 10% 6-4PP-CPD regions - with 64PP repair
T10_Diff_enh_w_repair<-T10_Diff_enh_XR64[which(T10_Diff_enh_XR64$Total_XR64_strand_mean_repair!=0),]
T10_Diff_enh_w_repair$Group<-"Top 10% 6-4PP-CPD - 64PP repair"

T10_Enh_Diff_XR64_mutR_table[3,1]<-nrow(T10_Diff_enh_w_repair) # All enhancers with 64PP repair in top 6-4PP-CPD regions
T10_Enh_Diff_XR64_mutR_table[3,2]<-nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with 64PP repair in top 64 with C>T mutations
T10_Enh_Diff_XR64_mutR_table[3,3]<-nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate==0),]) # Enh with 64PP repair in top 64 with NO C>T mutations
T10_Enh_Diff_XR64_mutR_table[3,4]<-(nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate!=0),])/nrow(T10_Diff_enh_w_repair))*100 # Percent of enh with 64PP repair in top 64 with C>T mutations
T10_Enh_Diff_XR64_mutR_table[3,5]<-(nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate==0),])/nrow(T10_Diff_enh_w_repair))*100 # Percent of enh with 64PP repair in top 64 with NO C>T mutations

# Remaining stats
T10_Enh_Diff_XR64_mutR_table[4,2]<-(T10_Enh_Diff_XR64_mutR_table[2,2]/T10_Enh_Diff_XR64_mutR_table[1,2])*100 # Percent mutated enhancers without 64PP repair
T10_Enh_Diff_XR64_mutR_table[5,2]<-(T10_Enh_Diff_XR64_mutR_table[3,2]/T10_Enh_Diff_XR64_mutR_table[1,2])*100 # Percent mutated enhancers with 64PP repair
T10_Enh_Diff_XR64_mutR_table[4,3]<-(T10_Enh_Diff_XR64_mutR_table[2,3]/T10_Enh_Diff_XR64_mutR_table[1,3])*100 # Percent enhancers not mutated without 64PP repair
T10_Enh_Diff_XR64_mutR_table[5,3]<-(T10_Enh_Diff_XR64_mutR_table[3,3]/T10_Enh_Diff_XR64_mutR_table[1,3])*100 # Percent enhancers not mutated with 64PP repair

write.csv(T10_Enh_Diff_XR64_mutR_table, file = paste0(out,"T10 6-4PP-CPD enhancers - XR64 repair and mutation stats - w ChrX.csv"))


# 6-4PP repair analysis - Bottom 10% 6-4PP-CPD -----------------------------

## TABLE with enhancer stats in bottom 6-4PP-CPD regions
B10_Enh_Diff_XR64_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(B10_Enh_Diff_XR64_mutR_table)<-c("All Enh in Bottom 6-4PP-CPD","Enhancers without 6-4PP repair","Enhancers with 6-4PP repair","Percent  without repair","Percent  with repair")
colnames(B10_Enh_Diff_XR64_mutR_table)<-c("All Enh in Bottom 6-4PP-CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Bottom 10%
B10_Diff_enh_XR64<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XR64,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% 6-4PP-CPD")
B10_Diff_enh_XR64$overlap<-NULL
B10_Diff_enh_XR64$percent_overlap<-NULL
B10_Diff_enh_XR64$Region<-NULL
B10_Diff_enh_XR64<-as.data.frame(B10_Diff_enh_XR64)

# Stats
B10_Enh_Diff_XR64_mutR_table[1,1]<-nrow(B10_Diff_enh_XR64) # All enhancers in bottom region
B10_Enh_Diff_XR64_mutR_table[1,2]<-nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
B10_Enh_Diff_XR64_mutR_table[1,3]<-nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
B10_Enh_Diff_XR64_mutR_table[1,4]<-(nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_C_to_T_rate!=0),])/nrow(B10_Diff_enh_XR64))*100 # Percent of all enhancers mutated
B10_Enh_Diff_XR64_mutR_table[1,5]<-(nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_C_to_T_rate==0),])/nrow(B10_Diff_enh_XR64))*100 # Percent of all enhancers NOT mutations

B10_Enh_Diff_XR64_mutR_table[4,1]<-(nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_XR64_strand_mean_repair==0),])/nrow(B10_Diff_enh_XR64))*100 # Percent enhancers without repair
B10_Enh_Diff_XR64_mutR_table[5,1]<-(nrow(B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_XR64_strand_mean_repair!=0),])/nrow(B10_Diff_enh_XR64))*100 # Percent enhancers with repair

# Enhancers in bottom 10% 6-4PP-CPD regions - no 64PP repair
B10_Diff_enh_no_repair<-B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_XR64_strand_mean_repair==0),]
B10_Diff_enh_no_repair$Group<-"Bottom 10% 6-4PP-CPD - no repair"

B10_Enh_Diff_XR64_mutR_table[2,1]<-nrow(B10_Diff_enh_no_repair) # All enhancers with no repair in bottom 6-4PP-CPD regions
B10_Enh_Diff_XR64_mutR_table[2,2]<-nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XR64_mutR_table[2,3]<-nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in bottom 6-4PP-CPD with NO C>T mutations
B10_Enh_Diff_XR64_mutR_table[2,4]<-(nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate!=0),])/nrow(B10_Diff_enh_no_repair))*100 # Percent of enh with no repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XR64_mutR_table[2,5]<-(nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate==0),])/nrow(B10_Diff_enh_no_repair))*100 # Percent of enh with no repair in bottom 6-4PP-CPD with NO C>T mutations

# Enhancers in bottom 10% 6-4PP-CPD regions - with 64PP repair
B10_Diff_enh_w_repair<-B10_Diff_enh_XR64[which(B10_Diff_enh_XR64$Total_XR64_strand_mean_repair!=0),]
B10_Diff_enh_w_repair$Group<-"Bottom 10% 6-4PP-CPD - 64PP repair"

B10_Enh_Diff_XR64_mutR_table[3,1]<-nrow(B10_Diff_enh_w_repair) # All enhancers with 64PP repair in bottom 6-4PP-CPD regions
B10_Enh_Diff_XR64_mutR_table[3,2]<-nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with 64PP repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XR64_mutR_table[3,3]<-nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate==0),]) # Enh with 64PP repair in bottom 6-4PP-CPD with NO C>T mutations
B10_Enh_Diff_XR64_mutR_table[3,4]<-(nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate!=0),])/nrow(B10_Diff_enh_w_repair))*100 # Percent of enh with 64PP repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XR64_mutR_table[3,5]<-(nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate==0),])/nrow(B10_Diff_enh_w_repair))*100 # Percent of enh with 64PP repair in bottom 6-4PP-CPD with NO C>T mutations

# Remaining stats
B10_Enh_Diff_XR64_mutR_table[4,2]<-(B10_Enh_Diff_XR64_mutR_table[2,2]/B10_Enh_Diff_XR64_mutR_table[1,2])*100 # Percent mutated enhancers without 64PP repair
B10_Enh_Diff_XR64_mutR_table[5,2]<-(B10_Enh_Diff_XR64_mutR_table[3,2]/B10_Enh_Diff_XR64_mutR_table[1,2])*100 # Percent mutated enhancers with 64PP repair
B10_Enh_Diff_XR64_mutR_table[4,3]<-(B10_Enh_Diff_XR64_mutR_table[2,3]/B10_Enh_Diff_XR64_mutR_table[1,3])*100 # Percent enhancers not mutated without 64PP repair
B10_Enh_Diff_XR64_mutR_table[5,3]<-(B10_Enh_Diff_XR64_mutR_table[3,3]/B10_Enh_Diff_XR64_mutR_table[1,3])*100 # Percent enhancers not mutated with 64PP repair
B10_Enh_Diff_XR64_mutR_table
write.csv(B10_Enh_Diff_XR64_mutR_table, file = paste0(out,"B10 6-4PP-CPD enhancers - XR64 repair and mutation stats - w ChrX.csv"))


# CPD repair analysis - Top 10% 6-4PP-CPD ----------------------------------

## TABLE with enhancer stats in top 6-4PP-CPD regions
T10_Enh_Diff_XRcpd_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(T10_Enh_Diff_XRcpd_mutR_table)<-c("All Enh in Top 6-4PP-CPD","Enhancers without CPD repair","Enhancers with CPD repair","Percent  without repair","Percent  with repair")
colnames(T10_Enh_Diff_XRcpd_mutR_table)<-c("All Enh in Top 6-4PP-CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Top 10%
T10_Diff_Enh_XRcpd<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = Diff_T10,Overlap = 50,Region_name = "Top 10% 6-4PP-CPD")
T10_Diff_Enh_XRcpd$overlap<-NULL
T10_Diff_Enh_XRcpd$percent_overlap<-NULL
T10_Diff_Enh_XRcpd$Region<-NULL
T10_Diff_Enh_XRcpd<-as.data.frame(T10_Diff_Enh_XRcpd)

# Stats
T10_Enh_Diff_XRcpd_mutR_table[1,1]<-nrow(T10_Diff_Enh_XRcpd) # All enhancers in top region
T10_Enh_Diff_XRcpd_mutR_table[1,2]<-nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
T10_Enh_Diff_XRcpd_mutR_table[1,3]<-nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
T10_Enh_Diff_XRcpd_mutR_table[1,4]<-(nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_C_to_T_rate!=0),])/nrow(T10_Diff_Enh_XRcpd))*100 # Percent of all enhancers mutated
T10_Enh_Diff_XRcpd_mutR_table[1,5]<-(nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_C_to_T_rate==0),])/nrow(T10_Diff_Enh_XRcpd))*100 # Percent of all enhancers NOT mutations

T10_Enh_Diff_XRcpd_mutR_table[4,1]<-(nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),])/nrow(T10_Diff_Enh_XRcpd))*100 # Percent enhancers without repair
T10_Enh_Diff_XRcpd_mutR_table[5,1]<-(nrow(T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),])/nrow(T10_Diff_Enh_XRcpd))*100 # Percent enhancers with repair

# Enhancers in top 10% 6-4PP-CPD regions - no CPD repair
T10_Diff_enh_no_repair<-T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),]
T10_Diff_enh_no_repair$Group<-"Top 10% 6-4PP-CPD - no repair"

T10_Enh_Diff_XRcpd_mutR_table[2,1]<-nrow(T10_Diff_enh_no_repair) # All enhancers with no repair in top 6-4PP-CPD regions
T10_Enh_Diff_XRcpd_mutR_table[2,2]<-nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[2,3]<-nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in top 6-4PP-CPD with NO C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[2,4]<-(nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate!=0),])/nrow(T10_Diff_enh_no_repair))*100 # Percent of enh with no repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[2,5]<-(nrow(T10_Diff_enh_no_repair[which(T10_Diff_enh_no_repair$Total_C_to_T_rate==0),])/nrow(T10_Diff_enh_no_repair))*100 # Percent of enh with no repair in top 6-4PP-CPD with NO C>T mutations

# Enhancers in top 10% 6-4PP-CPD regions - with CPD repair
T10_Diff_enh_w_repair<-T10_Diff_Enh_XRcpd[which(T10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),]
T10_Diff_enh_w_repair$Group<-"Top 10% 6-4PP-CPD - CPD repair"

T10_Enh_Diff_XRcpd_mutR_table[3,1]<-nrow(T10_Diff_enh_w_repair) # All enhancers with CPD repair in top 6-4PP-CPD regions
T10_Enh_Diff_XRcpd_mutR_table[3,2]<-nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with CPD repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[3,3]<-nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate==0),]) # Enh with CPD repair in top 6-4PP-CPD with NO C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[3,4]<-(nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate!=0),])/nrow(T10_Diff_enh_w_repair))*100 # Percent of enh with CPD repair in top 6-4PP-CPD with C>T mutations
T10_Enh_Diff_XRcpd_mutR_table[3,5]<-(nrow(T10_Diff_enh_w_repair[which(T10_Diff_enh_w_repair$Total_C_to_T_rate==0),])/nrow(T10_Diff_enh_w_repair))*100 # Percent of enh with CPD repair in top 6-4PP-CPD with NO C>T mutations

# Remaining stats
T10_Enh_Diff_XRcpd_mutR_table[4,2]<-(T10_Enh_Diff_XRcpd_mutR_table[2,2]/T10_Enh_Diff_XRcpd_mutR_table[1,2])*100 # Percent mutated enhancers without CPD repair
T10_Enh_Diff_XRcpd_mutR_table[5,2]<-(T10_Enh_Diff_XRcpd_mutR_table[3,2]/T10_Enh_Diff_XRcpd_mutR_table[1,2])*100 # Percent mutated enhancers with CPD repair
T10_Enh_Diff_XRcpd_mutR_table[4,3]<-(T10_Enh_Diff_XRcpd_mutR_table[2,3]/T10_Enh_Diff_XRcpd_mutR_table[1,3])*100 # Percent enhancers not mutated without CPD repair
T10_Enh_Diff_XRcpd_mutR_table[5,3]<-(T10_Enh_Diff_XRcpd_mutR_table[3,3]/T10_Enh_Diff_XRcpd_mutR_table[1,3])*100 # Percent enhancers not mutated with CPD repair
T10_Enh_Diff_XRcpd_mutR_table
write.csv(T10_Enh_Diff_XRcpd_mutR_table, file = paste0(out,"T10 6-4PP-CPD enhancers - XRcpd repair and mutation stats - w ChrX.csv"))


# CPD repair analysis - Bottom 10% 6-4PP-CPD -------------------------------

## TABLE with enhancer stats in bottom 6-4PP-CPD regions
B10_Enh_Diff_XRcpd_mutR_table<-as.data.frame( matrix(nrow = 5,ncol = 5))
rownames(B10_Enh_Diff_XRcpd_mutR_table)<-c("All Enh in Bottom 6-4PP-CPD","Enhancers without CPD repair","Enhancers with CPD repair","Percent  without repair","Percent  with repair")
colnames(B10_Enh_Diff_XRcpd_mutR_table)<-c("All Enh in Bottom 6-4PP-CPD","Enhancers with mutation rate > 0","Ehnancers with no mutation rate","Percent mutated", "Percent not mutated")

## Bottom 10%
B10_Diff_Enh_XRcpd<-Gene_Ovrlp_fun(Query = makeGRangesFromDataFrame(Enh_XRcpd,keep.extra.columns = T),Subject = Diff_B10,Overlap = 50,Region_name = "Bottom 10% 6-4PP-CPD")
B10_Diff_Enh_XRcpd$overlap<-NULL
B10_Diff_Enh_XRcpd$percent_overlap<-NULL
B10_Diff_Enh_XRcpd$Region<-NULL
B10_Diff_Enh_XRcpd<-as.data.frame(B10_Diff_Enh_XRcpd)

# Stats
B10_Enh_Diff_XRcpd_mutR_table[1,1]<-nrow(B10_Diff_Enh_XRcpd) # All enhancers in bottom region
B10_Enh_Diff_XRcpd_mutR_table[1,2]<-nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_C_to_T_rate!=0),]) # All enhancers with C>T mutation rate
B10_Enh_Diff_XRcpd_mutR_table[1,3]<-nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_C_to_T_rate==0),]) # All enhancers with NO C>T mutation rate
B10_Enh_Diff_XRcpd_mutR_table[1,4]<-(nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_C_to_T_rate!=0),])/nrow(B10_Diff_Enh_XRcpd))*100 # Percent of all enhancers mutated
B10_Enh_Diff_XRcpd_mutR_table[1,5]<-(nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_C_to_T_rate==0),])/nrow(B10_Diff_Enh_XRcpd))*100 # Percent of all enhancers NOT mutations

B10_Enh_Diff_XRcpd_mutR_table[4,1]<-(nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),])/nrow(B10_Diff_Enh_XRcpd))*100 # Percent enhancers without repair
B10_Enh_Diff_XRcpd_mutR_table[5,1]<-(nrow(B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),])/nrow(B10_Diff_Enh_XRcpd))*100 # Percent enhancers with repair

# Enhancers in bottom 10% 6-4PP-CPD regions - no CPD repair
B10_Diff_enh_no_repair<-B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair==0),]
B10_Diff_enh_no_repair$Group<-"Bottom 10% 6-4PP-CPD - no repair"

B10_Enh_Diff_XRcpd_mutR_table[2,1]<-nrow(B10_Diff_enh_no_repair) # All enhancers with no repair in bottom 6-4PP-CPD regions
B10_Enh_Diff_XRcpd_mutR_table[2,2]<-nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate!=0),]) # Enh with no repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[2,3]<-nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate==0),]) # Enh with no repair in bottom 6-4PP-CPD with NO C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[2,4]<-(nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate!=0),])/nrow(B10_Diff_enh_no_repair))*100 # Percent of enh with no repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[2,5]<-(nrow(B10_Diff_enh_no_repair[which(B10_Diff_enh_no_repair$Total_C_to_T_rate==0),])/nrow(B10_Diff_enh_no_repair))*100 # Percent of enh with no repair in bottom 6-4PP-CPD with NO C>T mutations

# Enhancers in bottom 10% 6-4PP-CPD regions - with CPD repair
B10_Diff_enh_w_repair<-B10_Diff_Enh_XRcpd[which(B10_Diff_Enh_XRcpd$Total_XRCPD_strand_mean_repair!=0),]
B10_Diff_enh_w_repair$Group<-"Bottom 10% 6-4PP-CPD - CPD repair"

B10_Enh_Diff_XRcpd_mutR_table[3,1]<-nrow(B10_Diff_enh_w_repair) # All enhancers with CPD repair in bottom 6-4PP-CPD regions
B10_Enh_Diff_XRcpd_mutR_table[3,2]<-nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate!=0),]) # Enh with CPD repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[3,3]<-nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate==0),]) # Enh with CPD repair in bottom 6-4PP-CPD with NO C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[3,4]<-(nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate!=0),])/nrow(B10_Diff_enh_w_repair))*100 # Percent of enh with CPD repair in bottom 6-4PP-CPD with C>T mutations
B10_Enh_Diff_XRcpd_mutR_table[3,5]<-(nrow(B10_Diff_enh_w_repair[which(B10_Diff_enh_w_repair$Total_C_to_T_rate==0),])/nrow(B10_Diff_enh_w_repair))*100 # Percent of enh with CPD repair in bottom 6-4PP-CPD with NO C>T mutations

# Remaining stats
B10_Enh_Diff_XRcpd_mutR_table[4,2]<-(B10_Enh_Diff_XRcpd_mutR_table[2,2]/B10_Enh_Diff_XRcpd_mutR_table[1,2])*100 # Percent mutated enhancers without CPD repair
B10_Enh_Diff_XRcpd_mutR_table[5,2]<-(B10_Enh_Diff_XRcpd_mutR_table[3,2]/B10_Enh_Diff_XRcpd_mutR_table[1,2])*100 # Percent mutated enhancers with CPD repair
B10_Enh_Diff_XRcpd_mutR_table[4,3]<-(B10_Enh_Diff_XRcpd_mutR_table[2,3]/B10_Enh_Diff_XRcpd_mutR_table[1,3])*100 # Percent enhancers not mutated without CPD repair
B10_Enh_Diff_XRcpd_mutR_table[5,3]<-(B10_Enh_Diff_XRcpd_mutR_table[3,3]/B10_Enh_Diff_XRcpd_mutR_table[1,3])*100 # Percent enhancers not mutated with CPD repair
B10_Enh_Diff_XRcpd_mutR_table
write.csv(B10_Enh_Diff_XRcpd_mutR_table, file = paste0(out,"B10 6-4PP-CPD enhancers - XRcpd repair and mutation stats - w ChrX.csv"))


### END ###
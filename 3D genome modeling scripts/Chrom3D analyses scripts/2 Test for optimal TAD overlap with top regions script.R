library(GenomicRanges)


## Description: This R script is meant to identify and select TADs used in the 3D genome model for IMR90 that overlap with the top 100 kb regions of interest (top 10%). 
# Will test different thresholds of overlap to determine which threshold has the optimal amount of coverage of top regions within TADs


out<-"output_directory/"


# TAD genomic positions for IMR90
imrTADS<-read.table("file_path/IMR90_inter_intra_chr_w_LADs.gtrack", header = F)
colnames(imrTADS)<-c("seqid","start","end","id","radius","periphery","edges" )
imrTADS$start<-imrTADS$start+1

TADS<-makeGRangesFromDataFrame(imrTADS, keep.extra.columns = T)


## Basepair overlap thresholds to measure: 1 bp (smallest amount of overlap), 10%, 20%, 30%, 40%, 50% overlap. Six measurements in total

# Table containing stat info for TOP TADs
# Parameters: 
param<-c('total_TADs',
         'TAD_median_size_WG',
         'TADs_with_any_overlap',
         'TADS_with_10_perc_ovlp',
         'TADS_with_20_perc_ovlp',
         'TADS_with_30_perc_ovlp',
         'TADS_with_40_perc_ovlp',
         'TADS_with_50_perc_ovlp',
         
         'perc_of_TADs_with_any_overlap',
         'perc_of_TADs_with_10%_overlap',
         'perc_of_TADs_with_20%_overlap',
         'perc_of_TADs_with_30%_overlap',
         'perc_of_TADs_with_40%_overlap',
         'perc_of_TADs_with_50%_overlap',
         
         'perc_coverage_of_Top_region_with_any_overlap',
         'perc_coverage_of_Top_region_with_10%_overlap',
         'perc_coverage_of_Top_region_with_20%_overlap',
         'perc_coverage_of_Top_region_with_30%_overlap',
         'perc_coverage_of_Top_region_with_40%_overlap',
         'perc_coverage_of_Top_region_with_50%_overlap',
         
         'Non_unique_TADs.any_ovrlp', # new
         'Non_unique_TADs.10%_ovrlp',
         'Non_unique_TADs.20%_ovrlp',
         'Non_unique_TADs.30%_ovrlp',
         'Non_unique_TADs.40%_ovrlp',
         'Non_unique_TADs.50%_ovrlp'
         
)

TAD_ovlp_stats<-as.data.frame(matrix(nrow = length(param), ncol = 6))
rownames(TAD_ovlp_stats)<-param
colnames(TAD_ovlp_stats)<-c('Top_10_64PP','Top_10_CPD','Top_10_RankDiff','Bottom_10_64PP','Bottom_10_CPD','Bottom_10_RankDiff' )

TAD_ovlp_stats[1,]<-nrow(imrTADS)
TAD_ovlp_stats[2,]<-median(TADS@ranges@width)


# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff, keep.extra.columns = T)

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))

# get total size of top regions
Top_64_bp<-sum(IP64_T10@ranges@width)
Bottom_64_bp<-sum(IP64_B10@ranges@width)

Top_CPD_bp<-sum(CPD_T10@ranges@width)
Bottom_CPD_bp<-sum(CPD_B10@ranges@width)

Top_diff_bp<-sum(Diff_T10@ranges@width)
Bottom_diff_bp<-sum(Diff_B10@ranges@width)

######
## Measure TAD overlaps using new function

Top_regions<-GRangesList(IP64_T10, CPD_T10, Diff_T10, IP64_B10, CPD_B10, Diff_B10)
Region_bp<-list(Top_64_bp,Top_CPD_bp, Top_diff_bp, Bottom_64_bp,Bottom_CPD_bp, Bottom_diff_bp)


# Function to measure variable TAD overlap
TAD_ovrlp_fun<-function(TAD_regions, Select_regions, Region_bp) {
  
  for (i in 1:length(Select_regions)) {
    Regions<-unlist(Select_regions[i])
    Top_bp<-unlist(Region_bp[i])
    
    TAD_ovrlps<-TAD_regions
    hits <- findOverlaps(TAD_regions,Regions)
    gr.over<- pintersect(TAD_regions[queryHits(hits)],Regions[subjectHits(hits)])
    gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
      sum(width(x)))
    TAD_ovrlps$overlap<-0
    TAD_ovrlps$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
    TAD_ovrlps<-data.frame(TAD_ovrlps)
    TAD_ovrlps$percent_overlap<-(TAD_ovrlps$overlap/TAD_ovrlps$width) *100
    
    # Measurement: Amount of TADs with N overlaps
    TAD_ovlp_stats[3,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$overlap >=1),])
    TAD_ovlp_stats[4,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=10),])
    TAD_ovlp_stats[5,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=20),])
    TAD_ovlp_stats[6,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=30),])
    TAD_ovlp_stats[7,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=40),])
    TAD_ovlp_stats[8,i]<-nrow(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=50),])
    
    ## Measurement: Percentage of TADs with overlaps with top regions ( rows 9:14)
    TAD_ovlp_stats[9,i]<-(TAD_ovlp_stats[3,i]/TAD_ovlp_stats[3,i])*100
    TAD_ovlp_stats[10,i]<-(TAD_ovlp_stats[4,i]/TAD_ovlp_stats[3,i])*100
    TAD_ovlp_stats[11,i]<-(TAD_ovlp_stats[5,i]/TAD_ovlp_stats[3,i])*100
    TAD_ovlp_stats[12,i]<-(TAD_ovlp_stats[6,i]/TAD_ovlp_stats[3,i])*100
    TAD_ovlp_stats[13,i]<-(TAD_ovlp_stats[7,i]/TAD_ovlp_stats[3,i])*100
    TAD_ovlp_stats[14,i]<-(TAD_ovlp_stats[8,i]/TAD_ovlp_stats[3,i])*100
    
    TAD_ovlp_stats[15,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$overlap >=1),10])/Top_bp)*100
    TAD_ovlp_stats[16,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=10),10])/Top_bp)*100
    TAD_ovlp_stats[17,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=20),10])/Top_bp)*100
    TAD_ovlp_stats[18,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=30),10])/Top_bp)*100
    TAD_ovlp_stats[19,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=40),10])/Top_bp)*100
    TAD_ovlp_stats[20,i]<-(sum(TAD_ovrlps[which(TAD_ovrlps$percent_overlap >=50),10])/Top_bp)*100
    
    
  }
  
    # Measure if TADs appear twice in both top and bottom regions depending on threshold used
  
  for (i in 1:3) { 
    Region_1<-unlist(Select_regions[i])
    Region_2<-unlist(Select_regions[i+3])
    
    TAD_ovrlps_1<-TAD_regions
    hits <- findOverlaps(TAD_regions,Region_1)
    gr.over<- pintersect(TAD_regions[queryHits(hits)],Region_1[subjectHits(hits)])
    gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
      sum(width(x)))
    TAD_ovrlps_1$overlap<-0
    TAD_ovrlps_1$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
    TAD_ovrlps_1<-data.frame(TAD_ovrlps_1)
    TAD_ovrlps_1$percent_overlap<-(TAD_ovrlps_1$overlap/TAD_ovrlps_1$width) *100
    
    TAD_ovrlps_2<-TAD_regions
    hits <- findOverlaps(TAD_regions,Region_2)
    gr.over<- pintersect(TAD_regions[queryHits(hits)],Region_2[subjectHits(hits)])
    gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
      sum(width(x)))
    TAD_ovrlps_2$overlap<-0
    TAD_ovrlps_2$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
    TAD_ovrlps_2<-data.frame(TAD_ovrlps_2)
    TAD_ovrlps_2$percent_overlap<-(TAD_ovrlps_2$overlap/TAD_ovrlps_2$width) *100
    
    # any overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$overlap >1),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$overlap >1),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[21,i]<-NROW(df1)
    
    # 10% overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$percent_overlap >=10),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$percent_overlap >=10),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[22,i]<-NROW(df1)
    
    # 20% overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$percent_overlap >=20),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$percent_overlap >=20),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[23,i]<-NROW(df1)
    
    # 30% overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$percent_overlap >=30),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$percent_overlap >=30),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[24,i]<-NROW(df1)
    
    # 40% overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$percent_overlap >=40),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$percent_overlap >=40),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[25,i]<-NROW(df1)
    
    # 50% overlap
    df1<-TAD_ovrlps_1[which(TAD_ovrlps_1$percent_overlap >=50),]
    df1<-makeGRangesFromDataFrame(df1)
    df2<-TAD_ovrlps_2[which(TAD_ovrlps_2$percent_overlap >=50),]
    df2<-makeGRangesFromDataFrame(df2)
    df1<-subsetByOverlaps(x = df1,ranges = df2)
    TAD_ovlp_stats[26,i]<-NROW(df1)
    
  }
  
  return(TAD_ovlp_stats)
}

TAD_ovlps_stats_TEST<-TAD_ovrlp_fun(TAD_regions = TADS, Select_regions = Top_regions, Region_bp = Region_bp)
write.csv(TAD_ovlps_stats_TEST, file = paste0(out, "TAD_stats_with_variable_overlap_100kb.csv"))



## END ##
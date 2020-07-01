library(GenomicRanges)


## Description: This R script is meant to select TADs that contain at least >10% overlap with 100kb regions that will be highlighted in the 3D genome (previously determined by measuring the overlap with various thresholds)


out<-"output_directory/"


# TAD genomic positions for IMR90
imrTADS<-read.table("file_path/IMR90_inter_intra_chr_w_LADs.gtrack", header = F)
colnames(imrTADS)<-c("seqid","start","end","id","radius","periphery","edges" )
imrTADS$start<-imrTADS$start+1
TADS<-makeGRangesFromDataFrame(imrTADS, keep.extra.columns = T)

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


######
## Top 10% regions 

## Calculate number of TADs overlaps. Only select TADs with >=10% overlap and generate TAD ID list

# function to subset overlaps between top regions and TAD coordinates
TAD_ovrlps<-function(TAD_regions, Select_regions ) {
  
  TAD_ovrlps<-TAD_regions
  hits <- findOverlaps(TAD_regions,Select_regions)
  gr.over<- pintersect(TAD_regions[queryHits(hits)],Select_regions[subjectHits(hits)])
  gr.counts <- tapply(gr.over,queryHits(hits),FUN=function(x)
    sum(width(x)))
  TAD_ovrlps$overlap<-0
  TAD_ovrlps$overlap[as.numeric(names(gr.counts))]<- unname(gr.counts)
  TAD_ovrlps<-data.frame(TAD_ovrlps)
  TAD_ovrlps$percent_overlap<-(TAD_ovrlps$overlap/TAD_ovrlps$width) *100
  TAD_ovrlps<-makeGRangesFromDataFrame(TAD_ovrlps,keep.extra.columns = T)
  return(TAD_ovrlps)
}

# Function to generate TAD ID list for both A & B chromosomes
TAD_ID<-function(TAD_regions ) {
  
  TAD_IDs<-TAD_regions$id
  TAD_IDs<-as.character(TAD_IDs)
  
  listA<-c()
  for (i in 1:length(TAD_IDs)) {
    list<-unlist(strsplit(TAD_IDs[i], ":"))
    list<-paste0(list[1], "_A:", list[2])
    listA<-append(listA, list)
  }
  listB<-c()
  for (i in 1:length(TAD_IDs)) {
    list<-unlist(strsplit(TAD_IDs[i], ":"))
    list<-paste0(list[1], "_B:", list[2])
    listB<-append(listB, list)
  }
  listAB<-append(listA,listB)
  listAB<-as.character(listAB)
  return(listAB)
}


## Top 10% 64PP - TAD overlaps (10%overlap) 
TADs_T64<-TAD_ovrlps(TAD_regions = TADS,Select_regions = IP64_T10)
TADs_T64<-subset(TADs_T64, percent_overlap >=10 )

## Bottom 10% 64PP - TAD overlaps (10%overlap) 
TADs_B64<-TAD_ovrlps(TAD_regions = TADS,Select_regions = IP64_B10)
TADs_B64<-subset(TADs_B64, percent_overlap >=10 )

## Find and exclude overlapping ranges between top and bottom
TADs_T64<-subsetByOverlaps(x = TADs_T64,ranges = TADs_B64,invert = T)
TADs_B64<-subsetByOverlaps(x = TADs_B64,ranges = TADs_T64,invert = T)

IP64T_listAB<-TAD_ID(TAD_regions = TADs_T64)
write.table(IP64T_listAB, file = paste0(out, "T10_64PP.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)

IP64B_listAB<-TAD_ID(TAD_regions = TADs_B64)
write.table(IP64B_listAB, file = paste0(out, "B10_64PP.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)


## Top 10% CPD - TAD overlaps (10%overlap) 
TADs_T_CPD<-TAD_ovrlps(TAD_regions = TADS,Select_regions = CPD_T10)
TADs_T_CPD<-subset(TADs_T_CPD,percent_overlap>=10)

## Bottom 10% CPD - TAD overlaps (10%overlap) 
TADs_B_CPD<-TAD_ovrlps(TAD_regions = TADS,Select_regions = CPD_B10)
TADs_B_CPD<-subset(TADs_B_CPD,percent_overlap>=10)

## Find and exclude overlapping ranges between top and bottom
TADs_T_CPD<-subsetByOverlaps(x = TADs_T_CPD,ranges = TADs_B_CPD,invert = T)
TADs_B_CPD<-subsetByOverlaps(x = TADs_B_CPD,ranges = TADs_T_CPD,invert = T)

CPD_T_listAB<-TAD_ID(TAD_regions = TADs_T_CPD)
write.table(CPD_T_listAB, file = paste0(out, "T10_CPD.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)

CPD_B_listAB<-TAD_ID(TAD_regions = TADs_B_CPD)
write.table(CPD_B_listAB, file = paste0(out, "B10_CPD.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)


## Top 10% RankDiff - TAD overlaps (10%overlap) 
TADs_T_Diff<-TAD_ovrlps(TAD_regions = TADS,Select_regions = Diff_T10)
TADs_T_Diff<-subset(TADs_T_Diff,percent_overlap>=10)

## Bottom 10% RankDiff - TAD overlaps (10%overlap) 
TADs_B_Diff<-TAD_ovrlps(TAD_regions = TADS,Select_regions = Diff_B10)
TADs_B_Diff<-subset(TADs_B_Diff,percent_overlap>=10)

## Find and exclude overlapping ranges between top and bottom
TADs_T_Diff<-subsetByOverlaps(x = TADs_T_Diff,ranges = TADs_B_Diff,invert = T)
TADs_B_Diff<-subsetByOverlaps(x = TADs_B_Diff,ranges = TADs_T_Diff,invert = T)

Diff_T_listAB<-TAD_ID(TAD_regions = TADs_T_Diff)
write.table(Diff_T_listAB, file = paste0(out, "T10_RankDiff.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)

Diff_B_listAB<-TAD_ID(TAD_regions = TADs_B_Diff)
write.table(Diff_B_listAB, file = paste0(out, "B10_RankDiff.100kb_10per_TAD_ovrlp.ids"), sep = "", col.names = F, row.names = F, quote = F)


## END ##
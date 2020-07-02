library(GenomicRanges)
library(rtracklayer)
library(annotatr)


## Description: Enrichment analysis of super enhancer elements for IMR90 (dbSUPER database) within top and
# bottom 10% regions


out<-"output_directory/"

# enrichment function
load(file = "path_to_enrichment_function/Enrich_function.RData")


######
## Load the UV lesion dataset then calculate top 10% regions of interest

# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")

## set dataframe with ranges for whole genome
Genome<-Ranked_diff
Genome<-makeGRangesFromDataFrame(Genome)
WG_size<-sum(Ranked_diff$width)

Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff, keep.extra.columns = T)

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


######
## Load super enhancers dataset

# Super enhancers
SuEnh<-import.bed("file_path/IMR90_dbSUPER_super_enhancer_regions.bed")
SuEnh$annot.type<-"Super_enhancer"
SuEnh<-as.data.frame(SuEnh)
SuEnh$score<-NULL

WG_ALL_Annots<-SuEnh
Annots<-"Super_enhancer"

######
## Measure enrichment of super enhancers

## Top 10%
# Top - 64PP
IP64_T10_Enrich<-Enrich(target = WG_ALL_Annots, gr = IP64_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(IP64_T10_Enrich, file = paste0(out,"IP64_T10_SUPER_Enh_enrich_100kb.RData"))

# Top - CPD
CPD_T10_Enrch<-Enrich(target = WG_ALL_Annots, gr = CPD_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(CPD_T10_Enrch, file = paste0(out,"CPD_T10_SUPER_Enh_enrich_100kb.RData"))

# Top - RankDiff
Diff_T10_Enrch<-Enrich(target = WG_ALL_Annots, gr = Diff_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(Diff_T10_Enrch, file = paste0(out,"RankDiff_T10_SUPER_Enh_enrich_100kb.RData"))


## Bottom 10%
# Bottom - 64PP
IP64_B10_Enrich<-Enrich(target = WG_ALL_Annots, gr = IP64_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(IP64_B10_Enrich, file = paste0(out,"IP64_B10_SUPER_Enh_enrich_100kb.RData"))

# Bottom - CPD
CPD_B10_Enrch<-Enrich(target = WG_ALL_Annots, gr = CPD_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(CPD_B10_Enrch, file = paste0(out,"CPD_B10_SUPER_Enh_enrich_100kb.RData"))

# Bottom -RankDiff
Diff_B10_Enrch<-Enrich(target = WG_ALL_Annots, gr = Diff_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(Diff_B10_Enrch, file = paste0(out,"RankDiff_B10_SUPER_Enh_enrich_100kb.RData"))


## END ##
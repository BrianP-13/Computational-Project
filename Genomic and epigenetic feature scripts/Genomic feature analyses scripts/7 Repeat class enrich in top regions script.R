library(GenomicRanges)

## Description: Enrichment analysis of repetitive elements within top 10% 100kb regions of RankDiff


# Strategy: 
# Step 1: Load RankDiff file and subset to top 10% regions
# Step 2: Load RepeatMasker database (already formated) along with the list containing repeat classes and subFamilies

# repClass:
# Step 3: Create list to contain the enrichment results
# Step 4: Conduct repeat enrichment analysis


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
Ranked_diff

## Top and bottom 10% regions
IP64_T10<-subset(Ranked_diff,ranked_64PP > quantile(Ranked_diff$ranked_64PP,0.9))
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))

CPD_T10<-subset(Ranked_diff,ranked_CPD > quantile(Ranked_diff$ranked_CPD,0.9))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))

Diff_T10<-subset(Ranked_diff,Difference > quantile(Ranked_diff$Difference,0.9))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


######
## Load RepeatMasker database (already formated) along with the list containing repeat classes and subFamilies

# Import RepeatMasker database of repeat elements
load(file = "file_path/modified.RMSK.hg19.RData")
load(file = "file_path/Repeat Classes with subFamily list.RData")
str(Repeat_Class_list)


## Enrichment of repeat Class types

# New list of repeat Class annotations to use
rClass_annots<-names(Repeat_Class_list)

Regions_of_int<-GRangesList(IP64_T10,CPD_T10,Diff_T10,IP64_B10,CPD_B10,Diff_B10)
Region_names<-c('IP64_T10','CPD_T10','Diff_T10','IP64_B10','CPD_B10','Diff_B10')

# start of loop
rClass.Enrich_list<-list()

for (i in 1:6) {
  Region_T10<-unlist(Regions_of_int[i])
  Region_Enrich<-Enrich(target = rmsk_df,gr = Region_T10,annot.type = "repClass",Annots = rClass_annots,Population_size = WG_size)
  Enrich_df<-paste0('RepClass_Enrch_', Region_names[i]) 
  rClass.Enrich_list[[Enrich_df]]<-Region_Enrich
  
  print(paste0("Done with ",Region_names[i]))
}
rClass.Enrich_list
str(rClass.Enrich_list)

save(rClass.Enrich_list, file = paste0(out, "RepClass_Enrich_Top10perc.100kb.regions.RData"))



## END ##
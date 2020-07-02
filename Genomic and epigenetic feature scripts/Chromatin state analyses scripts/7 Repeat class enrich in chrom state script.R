library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Description: Will measure the enrichment of specific repeat elements, Class subtype, within each of the 15 chromatin states

## Strategy:
# Step 1: Load chromatin state BED file and repeat masker database 
# Step 2: Load previously formated repeatmasker dataframe and repeat class list

# Step 3: Make list of repeat Class types within repeat Class to use for enrichment function
# Step 4: Create list to contain the enrichment results for each chromatin state
# Step 5: Create for loop to go through each chromatin state and measure repeat enrichment. Enrichment results for each State will be saved in a list class object.


out<-"output_directory/"


# enrichment function
load(file = "path_to_enrichment_function/Enrich_function.RData")

# Load IMR90 15 core chromatin state dataset
core15<-file.path("file_path/E017_15_coreMarks_dense.bed.bgz")
Chrom15<-import.bed(core15)
mcols(Chrom15)$score<-NULL
mcols(Chrom15)$itemRgb<-NULL
mcols(Chrom15)$thick<-NULL
df_core15<-df_core15[which(df_core15@seqnames != "chrY"),]
df_core15<-df_core15[which(df_core15@seqnames != "chrM"),]


Genome<-df_core15
WG_size<-sum(Genome@ranges@width)

core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_names2<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF_Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")


# Import RepeatMasker database of repeat elements
load(file = "file_path/modified.RMSK.hg19.RData")
load(file = "file_path/Repeat Classes with subFamily list.RData")
str(Repeat_Class_list)

# List of repeat Class annotations to use
rClass_annots<-names(Repeat_Class_list)
rClass_annots

# start of loop
Enrich_list<-list()

for (i in 1:length(core15_names)) {
  State_df<-subset(df_core15, name == core15_names[i])
  State_Enrich<-Enrich(target = rmsk_df,gr = State_df,annot.type = "repClass",Annots = rClass_annots, Population_size = WG_size)
  
  Enrich_df<-paste0('RepClass_Enrch_', core15_names2[i]) 
  Enrich_list[[Enrich_df]]<-State_Enrich
  
  print(paste0("Done with ",core15_names2[i]))
}
Enrich_list
str(Enrich_list)

save(Enrich_list, file = paste0(out, "RepeatClass_15ChromState_Enrichment.RData"))


## END ##
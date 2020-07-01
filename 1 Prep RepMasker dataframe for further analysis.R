library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Description: Will modify and adjust the Repeat Masker dataframe for enrichment analyses


out<-"output_directory/"


# Import RepeatMasker database of repeat elements
repM<-file.path("path_to_RepeatMasker_dataset/hg19.database.rmsk.EXTRA.COLS.csv")
rmsk<-read.csv(repM, header = TRUE, sep = ",")
rmsk<-rmsk[,2:19]
rmsk<-makeGRangesFromDataFrame(rmsk, keep.extra.columns = TRUE)

rmsk<-keepSeqlevels(rmsk, value = standardChromosomes(rmsk), pruning.mode = "coarse")
rmsk<-sortSeqlevels(rmsk)
rmsk<-sort(rmsk)
rmsk_df<-as.data.frame(rmsk)


## repFamily types to remove, mostly ambiguous repeat annotations
Removal_List<-c("Unknown","DNA?", "hAT?", "Gypsy?", "TcMar?", "SINE?", "PiggyBac?", "Helitron?","ERVL?", "L1?", "Penelope?","Unknown?", "LTR?","Other")

for (i in 1:14) { 
  rmsk_df<-subset(rmsk_df, repFamily != Removal_List[i])
}

## Replace a couple repeat Family types with more specific names to better describe the repeat element
rmsk_df[,14]<-as.character(rmsk_df[,14])
rmsk_df[,13]<-as.character(rmsk_df[,13])

SINE_tRNA<-rmsk_df[which(rmsk_df$repClass=="SINE"),]
SINE_tRNA<-SINE_tRNA[which(SINE_tRNA$repFamily=="tRNA"),]

SINEs<-rmsk_df[which(rmsk_df$repClass=="SINE"),]
SINEs<-SINEs[which(SINEs$repFamily=="SINE"),]

rmsk_df[as.vector(rownames(SINE_tRNA)),14]<-"SINE2_tRNA"
rmsk_df[as.vector(rownames(SINEs)),14]<-"LFSINE_Vert"

# Rename the "DNA" class to "DNA TE"
DNA_Rep<-rmsk_df[which(rmsk_df$repClass=="DNA"),]
rmsk_df[as.vector(rownames(DNA_Rep)),13]<-"DNA TE"

# Change the Rolling Circle class to the broader DNA TE class
RC<-rmsk_df[which(rmsk_df$repClass=="RC"),]
rmsk_df[as.vector(rownames(RC)),13]<-"DNA TE"

# Change all RNA repeat Class to simply "RNA"
srpRNA<-rmsk_df[which(rmsk_df$repClass=="srpRNA"),]
snRNA<-rmsk_df[which(rmsk_df$repClass=="snRNA"),]
rRNA<-rmsk_df[which(rmsk_df$repClass=="rRNA"),]
scRNA<-rmsk_df[which(rmsk_df$repClass=="scRNA"),]
tRNA<-rmsk_df[which(rmsk_df$repClass=="tRNA"),]

rmsk_df[as.vector(rownames(srpRNA)),13]<-"RNA"
rmsk_df[as.vector(rownames(snRNA)),13]<-"RNA"
rmsk_df[as.vector(rownames(rRNA)),13]<-"RNA"
rmsk_df[as.vector(rownames(scRNA)),13]<-"RNA"
rmsk_df[as.vector(rownames(tRNA)),13]<-"RNA"


# Make list of repeat Family types within repeat Class
rClass_annots<-as.character(unique(rmsk_df$repClass))
Repeat_Class_list<-list()

for (i in 1:length(rClass_annots)) {
  df<-rmsk_df[which(rmsk_df$repClass == rClass_annots[i]),]
  RepFamily_types<-as.character(unique(df$repFamily))
  Repeat_Class_list[[rClass_annots[i]]]<-RepFamily_types
  
}

## Reorder the repeat class types in the list to include Non-LTR retrotransposons, then LTR retrotrans. then DNA TE, then RNA repeats, then basic DNA repeats
Repeat_Class_list<-Repeat_Class_list[c(2,4,3,5,8,7,1,6)]

# Reorder LINE elements
Repeat_Class_list$LINE<-Repeat_Class_list$LINE[c(3,1,2,4,5,6)]

# Reorder SINE elements
Repeat_Class_list$SINE<-Repeat_Class_list$SINE[c(2,1,5,3,4)]


save(Repeat_Class_list, file = paste0(out, "Repeat Classes with subFamily list.RData"))
save(rmsk_df, file = paste0(out,"modified.RMSK.hg19.RData"))


## END ##
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(annotatr)


## Description: Will measure the enrichment of genic features and enhancer elements within the 15 defined chromatin states

out<-"output_directory/"


# enrichment function
load(file = "path_to_enrichment_function/Enrich_function.RData")


# Load IMR90 15 core chromatin state dataset
core15<-file.path("file_path/E017_15_coreMarks_dense.bed.bgz")
df_core15<-import.bed(core15)
mcols(df_core15)$score<-NULL
mcols(df_core15)$itemRgb<-NULL
mcols(df_core15)$thick<-NULL
df_core15<-df_core15[which(df_core15@seqnames != "chrY"),]
df_core15<-df_core15[which(df_core15@seqnames != "chrM"),]
df_core15<-sortSeqlevels(df_core15)
df_core15<-sort(df_core15)

Genome<-granges(df_core15)
WG_size<-sum(Genome@ranges@width)

core15_names<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
core15_names2<-c("1_TssA","2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF_Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")


######
## Step 2: Select annotations and build annotated genomes

# annotation names for basicgenes
BGenes<-build_annotations(genome = 'hg19', annotations = "hg19_basicgenes")
interg<-build_annotations(genome = 'hg19', annotations = "hg19_genes_intergenic")

### Annotate whole genome

# annotate genic regions in whole genome
WG_GenicF<-annotate_regions(regions = Genome, annotations = BGenes,
                            minoverlap = 1L, ignore.strand = T, quiet = F)
WG_GenicF<-data.frame(WG_GenicF)
WG_GenicF<-WG_GenicF[,c(6:15)]
WG_GenicF<-WG_GenicF[which(WG_GenicF$annot.gene_id != 'NA'),] # remove regions with no gene ID info
WG_GenicF$annot.tx_id<-NULL
WG_GenicF$annot.gene_id<-NULL
WG_GenicF$annot.symbol<-NULL

# annotate intergenic regions in whole genome
WG_InterGenes<-annotate_regions(regions = Genome, annotations = interg, 
                                minoverlap = 1L, ignore.strand = T, quiet = F)
WG_InterGenes<-data.frame(WG_InterGenes)
WG_InterGenes<-WG_InterGenes[,c(6:15)]
WG_InterGenes$annot.tx_id<-NULL
WG_InterGenes$annot.gene_id<-NULL
WG_InterGenes$annot.symbol<-NULL


## Load IMR90 enhancer dataset
IMR90_Enh_pos<-import.bed(con = "file_path/IMR90_enhancer_list_EnhancerAtlas2.0.bed")
IMR90_Enh_pos<-as.data.frame(IMR90_Enh_pos)
IMR90_Enh_pos$start<-IMR90_Enh_pos$start-1
IMR90_Enh_pos<-makeGRangesFromDataFrame(IMR90_Enh_pos,keep.extra.columns = T)
IMR90_Enh_pos$name<-paste0("EnhAtlas-IMR90 enhancer ",seq(from=1,to=NROW(IMR90_Enh_pos)))
IMR90_Enh_pos

# Annotate whole genome
WG_enh<-IMR90_Enh_pos
WG_enh<-as.data.frame(WG_enh)
WG_enh$annot.type<-"Enhancer-EnhancerAtlas"

# Super enhancers
SuEnh<-import.bed("file_path/IMR90_dbSUPER_super_enhancer_regions.bed")
SuEnh$annot.type<-"Super_enhancer"
SuEnh<-as.data.frame(SuEnh)
SuEnh$score<-NULL

# Combine all ranges into one data frame
colnames(SuEnh)<-colnames(WG_GenicF)
colnames(WG_enh)<-colnames(WG_GenicF)

WG_ALL_Annots<-rbind(WG_GenicF, WG_InterGenes,WG_enh,SuEnh)
Annots<-unique(WG_ALL_Annots$annot.type)


## start of loop
Enrich_list<-list()

for (i in 1:length(core15_names)) {
  State_df<-subset(df_core15, name == core15_names[i])
  State_Enrich<-Enrich(target = WG_ALL_Annots,gr = State_df,annot.type = "annot.type",Annots = Annots, Population_size = WG_size)
  
  Enrich_df<-paste0('Genic.Enh_Enrch_', core15_names2[i]) 
  Enrich_list[[Enrich_df]]<-State_Enrich
  
  print(paste0("Done with ",core15_names2[i]))
}
str(Enrich_list)
Enrich_list

save(Enrich_list, file = paste0(out, "Genic.EnhAltas_15ChromState_Enrchment.RData"))


## END ##
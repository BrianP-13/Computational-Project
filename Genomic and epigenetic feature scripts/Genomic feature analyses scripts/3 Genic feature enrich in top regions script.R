library(GenomicRanges)
library(annotatr)


## Description: Enrichment analysis of genic features within top and bottom 10% regions of interest
# 100kb bin size for analysis


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
## Select annotations and build annotated genomes

# annotation names for basicgenes
Annots<-c("hg19_genes_promoters","hg19_genes_1to5kb", "hg19_genes_5UTRs", "hg19_genes_exons", "hg19_genes_introns", "hg19_genes_3UTRs")
IGR<-"hg19_genes_intergenic"
basicgenes = 'hg19_basicgenes'

BGenes<-build_annotations(genome = 'hg19', annotations = basicgenes)
interg<-build_annotations(genome = 'hg19', annotations = "hg19_genes_intergenic")

## Annotate whole genome
# annotate genic regions in whole genome
WG_GenicF<-annotate_regions(regions = Genome, annotations = BGenes,
                            minoverlap = 1L, ignore.strand = T, quiet = F)
WG_GenicF<-data.frame(WG_GenicF)
WG_GenicF<-WG_GenicF[,c(6:15)]
WG_GenicF<-WG_GenicF[which(WG_GenicF$annot.gene_id != 'NA'),] # remove regions with no gene ID info

# annotate intergenic regions in whole genome
WG_InterGenes<-annotate_regions(regions = Genome, annotations = interg, 
                                minoverlap = 1L, ignore.strand = T, quiet = F)
WG_InterGenes<-data.frame(WG_InterGenes)
WG_InterGenes<-WG_InterGenes[,c(6:15)]

# Combine all ranges into one data frame
WG_ALL_Annots<-rbind(WG_GenicF, WG_InterGenes)
Annots<-unique(WG_ALL_Annots$annot.type)


######
## Measure enrichment of genic features 

## Top 10%
# Top - 64PP
IP64_T10_Enrich<-Enrich(target = WG_ALL_Annots, gr = IP64_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(IP64_T10_Enrich, file = paste0(out,"IP64_T10_GenicF_enrich_100kb.RData"))

# Top - CPD
CPD_T10_Enrch<-Enrich(target = WG_ALL_Annots, gr = CPD_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(CPD_T10_Enrch, file = paste0(out,"CPD_T10_GenicF_enrich_100kb.RData"))

# Top - RankDiff
Diff_T10_Enrch<-Enrich(target = WG_ALL_Annots, gr = Diff_T10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(Diff_T10_Enrch, file = paste0(out,"RankDiff_T10_GenicF_enrich_100kb.RData"))


## Bottom 10%
# Bottom - 64PP
IP64_B10_Enrich<-Enrich(target = WG_ALL_Annots, gr = IP64_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(IP64_B10_Enrich, file = paste0(out,"IP64_B10_GenicF_enrich_100kb.RData"))

# Bottom - CPD
CPD_B10_Enrch<-Enrich(target = WG_ALL_Annots, gr = CPD_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(CPD_B10_Enrch, file = paste0(out,"CPD_B10_GenicF_enrich_100kb.RData"))

# Bottom -RankDiff
Diff_B10_Enrch<-Enrich(target = WG_ALL_Annots, gr = Diff_B10, annot.type = "annot.type", Annots = Annots,Population_size = WG_size)
save(Diff_B10_Enrch, file = paste0(out,"RankDiff_B10_GenicF_enrich_100kb.RData"))


## END ##
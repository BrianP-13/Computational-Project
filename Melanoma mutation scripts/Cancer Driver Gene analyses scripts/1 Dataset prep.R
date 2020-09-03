library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(biomaRt)


out<-"output_directory/"

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Function to get additional gene info
Gene_ID_fun<-function(Gene_ID_df) { 
  Gene_ID_df$hgnc_id<-NA
  Gene_ID_df$hgnc_symbol<-NA
  Gene_ID_df$entrezgene_id<-NA
  
  IDs<-Gene_ID_df$gene_id
  ID_aliases<-getBM(attributes = c("ensembl_gene_id","hgnc_id","hgnc_symbol","entrezgene_id"), 
                    filters = 'ensembl_gene_id', 
                    values = IDs, 
                    mart = ensembl)
  
  for (i in 1:length(Gene_ID_df$gene_id)) {
    ID<-Gene_ID_df$gene_id[i]
    new_IDs<-ID_aliases[ID_aliases$ensembl_gene_id %in% ID ,]
    
    if ( nrow(new_IDs)==1 ) {
      Gene_ID_df[i,which(colnames(Gene_ID_df)=="hgnc_id")]<-new_IDs[,2]
      Gene_ID_df[i,which(colnames(Gene_ID_df)=="hgnc_symbol")]<-new_IDs[,3]
      Gene_ID_df[i,which(colnames(Gene_ID_df)=="entrezgene_id")]<-new_IDs[,4]
    }
    
  }
  return(Gene_ID_df)
}

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


######
## PART 1: Prep new dataset to consolidate all UV lesion, repair and mutation measurements for protein coding genes

# load rank normalized UV lesion files with ranked difference column
load(file = "file_path/UV lesions rank normalized to norm dist 100kb.RData")
Ranked_diff<-makeGRangesFromDataFrame(Ranked_diff, keep.extra.columns = T)

## Top and bottom 10% regions
IP64_B10<-subset(Ranked_diff,ranked_64PP < quantile(Ranked_diff$ranked_64PP,0.1))
CPD_B10<-subset(Ranked_diff,ranked_CPD < quantile(Ranked_diff$ranked_CPD,0.1))
Diff_B10<-subset(Ranked_diff,Difference < quantile(Ranked_diff$Difference,0.1))


# Load melanoma cancer driver list with mutation frequency info
load(file = 'file_path/MELA_driver_list_with_mut_freq_info.RData')
# Load 64PP dataset with dmg and repair
load(file = "file_path/Gene_list_64PP_dmg_repair_mutR.RData")
# Load CPD dataset with dmg and repair
load(file = "file_path/Gene_list_CPD_dmg_repair_mutR.RData")
# Combine datasets together
UV_dmg_repair_mutR<-IP64_dmg_repair_mutR
mcols(UV_dmg_repair_mutR)<-cbind(mcols(UV_dmg_repair_mutR), mcols(CPD_dmg_repair_mutR)[c(3:5)])
mcols(UV_dmg_repair_mutR)<-mcols(UV_dmg_repair_mutR)[c(1:4,8:9,5,10,6,7)]
mcols(UV_dmg_repair_mutR)

UV_dmg_repair_mutR<-subset(UV_dmg_repair_mutR,seqnames!="chrM")
UV_dmg_repair_mutR<-subset(UV_dmg_repair_mutR,seqnames!="chrY")


###
## UV ratio measurements

## 64PP ratio
# Load 64PP dataset with dmg and repair
load(file = "file_path/Gene_list_64PP_dmg_repair_mutR.RData")
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,seqnames!="chrM")
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,seqnames!="chrY")
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,seqnames!="chrX")

# Next, remove all genes with zero repair and genes with Na/NaN UV damage values
IP64_dmg_repair_mutR<-as.data.frame(IP64_dmg_repair_mutR)
IP64_dmg_repair_mutR<-subset(IP64_dmg_repair_mutR,Total_XR64_strand_mean_repair !=0)

# Measure the ratio between 64PP FC and repair
IP64_dmg_repair_mutR$FC_repair_ratio<-IP64_dmg_repair_mutR$Mean_64PP_FC/IP64_dmg_repair_mutR$Total_XR64_strand_mean_repair
IP64_dmg_repair_mutR$Log2_UV_ratio<-log2(IP64_dmg_repair_mutR$FC_repair_ratio) # Log2 transform the UV lesion ratio
IP64_dmg_repair_mutR<-NaRV.omit(IP64_dmg_repair_mutR)
IP64_dmg_repair_mutR$Group<-"All genes"
head(IP64_dmg_repair_mutR)


## CPD ratio
# Load CPD dataset with dmg and repair
load(file = "file_path/Gene_list_CPD_dmg_repair_mutR.RData")
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,seqnames!="chrM")
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,seqnames!="chrY")
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,seqnames!="chrX")

# Next, remove all genes with zero repair and genes with Na/NaN UV damage values
CPD_dmg_repair_mutR<-as.data.frame(CPD_dmg_repair_mutR)
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,Total_XRCPD_strand_mean_repair !=0)

# Measure the ratio between CPD FC and repair
CPD_dmg_repair_mutR$FC_repair_ratio<-CPD_dmg_repair_mutR$Mean_CPD_FC/CPD_dmg_repair_mutR$Total_XRCPD_strand_mean_repair
CPD_dmg_repair_mutR$Log2_UV_ratio<-log2(CPD_dmg_repair_mutR$FC_repair_ratio) # Log2 transform the UV lesion ratio
CPD_dmg_repair_mutR<-NaRV.omit(CPD_dmg_repair_mutR)
CPD_dmg_repair_mutR$Group<-"All genes"
head(CPD_dmg_repair_mutR)



## Add UV ratio scores to master dataframe
UV_dmg_repair_mutR$UV64_ratio<-0
UV_dmg_repair_mutR$Log2_UV64_ratio<-0

UV_dmg_repair_mutR$UVcpd_ratio<-0
UV_dmg_repair_mutR$Log2_UVcpd_ratio<-0

UV_dmg_repair_mutR
NROW(UV_dmg_repair_mutR)



for (i in 1:NROW(UV_dmg_repair_mutR)){
  ID<-UV_dmg_repair_mutR$gene_id[i]
  
  if (nrow(IP64_dmg_repair_mutR[which(IP64_dmg_repair_mutR$gene_id==ID),])==1) {
    UV_dmg_repair_mutR$UV64_ratio[i]<-IP64_dmg_repair_mutR[which(IP64_dmg_repair_mutR$gene_id==ID),13]
    UV_dmg_repair_mutR$Log2_UV64_ratio[i]<-IP64_dmg_repair_mutR[which(IP64_dmg_repair_mutR$gene_id==ID),14]
  } else {
    UV_dmg_repair_mutR$UV64_ratio[i]<-NA
    UV_dmg_repair_mutR$Log2_UV64_ratio[i]<-NA
  }
  
  if (nrow(CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$gene_id==ID),])==1) {
    UV_dmg_repair_mutR$UVcpd_ratio[i]<-CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$gene_id==ID),13]
    UV_dmg_repair_mutR$Log2_UVcpd_ratio[i]<-CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$gene_id==ID),14]
  } else {
    UV_dmg_repair_mutR$UVcpd_ratio[i]<-NA
    UV_dmg_repair_mutR$Log2_UVcpd_ratio[i]<-NA
  }
  
  if (i==(NROW(UV_dmg_repair_mutR)/4)*1){
    print("25% done")
  }
  if (i==(NROW(UV_dmg_repair_mutR)/4)*2){
    print("50% done")
  }
  if (i==(NROW(UV_dmg_repair_mutR)/4)*3){
    print("75% done")
  }
  
}
UV_dmg_repair_mutR



UV_dmg_repair_mutR<-Gene_ID_fun(Gene_ID_df = as.data.frame(UV_dmg_repair_mutR))
head(UV_dmg_repair_mutR)

save(UV_dmg_repair_mutR,file = paste0(out, "Prot_coding_genes_w_repair_mutR_UVratio_info.RData"))




######
## PART 2: Prep new dataset to consolidate all UV lesion, repair and mutation measurements for enhancer regions

# Load enhancer datasets: 6-4PP XRseq, CPD XRseq and C>T mutation rate data 
load(file = 'file_path/Enh_Total_XR_and_mutR_dataset.RData')
Enh_XR_mutR<-makeGRangesFromDataFrame(Enh_XR_mutR,keep.extra.columns = T)
Enh_XR_mutR<-sortSeqlevels(Enh_XR_mutR)
Enh_XR_mutR<-sort(Enh_XR_mutR)
Enh_XR_mutR$Log10_mut_rate<-log10(Enh_XR_mutR$Total_C_to_T_rate)
Enh_XR_mutR


###
## 64PP ratio

# Load enhancer datasets: 6-4PP FC
load(file='file_path/IP64_FC_per_Enh.RData')
IP64_Enh

Enh_XR_mutR_UV64<-Enh_XR_mutR[,c(1:2,4:5)]
mcols(Enh_XR_mutR_UV64)<-cbind(mcols(Enh_XR_mutR_UV64), mcols(IP64_Enh)[2])
mcols(Enh_XR_mutR_UV64)

Enh_XR_mutR_UV64<-subset(Enh_XR_mutR_UV64,seqnames!="chrX")
Enh_XR_mutR_UV64

# Next, remove all enhancers with zero repair and enhancers with Na/NaN UV damage values
IP64_ratio_MutR<-as.data.frame(Enh_XR_mutR_UV64)
IP64_ratio_MutR<-subset(IP64_ratio_MutR,Total_XR64_strand_mean_repair !=0)
head(IP64_ratio_MutR)
hist(IP64_ratio_MutR$Total_XR64_strand_mean_repair)

# Measure the ratio between 64PP FC and repair
IP64_ratio_MutR$FC_repair_ratio<-IP64_ratio_MutR$Mean_64PP_FC/IP64_ratio_MutR$Total_XR64_strand_mean_repair
IP64_ratio_MutR$Log2_UV_ratio<-log2(IP64_ratio_MutR$FC_repair_ratio) # Log2 transform the UV lesion ratio
IP64_ratio_MutR<-NaRV.omit(IP64_ratio_MutR)
#IP64_ratio_MutR$Group<-"All enhancers"
head(IP64_ratio_MutR)

###
## CPD ratio
# Load enhancer datasets: 6-4PP FC
load(file='file_path/CPD_FC_per_Enh.RData')
CPD_Enh

Enh_XR_mutR_UVcpd<-Enh_XR_mutR[,c(1,3:5)]
mcols(Enh_XR_mutR_UVcpd)<-cbind(mcols(Enh_XR_mutR_UVcpd), mcols(CPD_Enh)[2])
mcols(Enh_XR_mutR_UVcpd)

Enh_XR_mutR_UVcpd<-subset(Enh_XR_mutR_UVcpd,seqnames!="chrX")
Enh_XR_mutR_UVcpd

# Next, remove all enhancers with zero repair and enhancers with Na/NaN UV damage values
CPD_dmg_repair_mutR<-as.data.frame(Enh_XR_mutR_UVcpd)
CPD_dmg_repair_mutR<-subset(CPD_dmg_repair_mutR,Total_XRCPD_strand_mean_repair !=0)

# Measure the ratio between CPD FC and repair
CPD_dmg_repair_mutR$FC_repair_ratio<-CPD_dmg_repair_mutR$Mean_CPD_FC/CPD_dmg_repair_mutR$Total_XRCPD_strand_mean_repair
CPD_dmg_repair_mutR$Log2_UV_ratio<-log2(CPD_dmg_repair_mutR$FC_repair_ratio) # Log2 transform the UV lesion ratio
CPD_dmg_repair_mutR<-NaRV.omit(CPD_dmg_repair_mutR)
#CPD_dmg_repair_mutR$Group<-"All enhancers"
head(CPD_dmg_repair_mutR)


### Create master table

# Combine datasets with repair, mutation and susc. together
UV_dmg_repair_mutR<-Enh_XR_mutR
mcols(UV_dmg_repair_mutR)<-cbind(mcols(UV_dmg_repair_mutR),mcols(IP64_Enh)[2],mcols(CPD_Enh)[2])
UV_dmg_repair_mutR

# measure log2(UV lesion FC)
UV_dmg_repair_mutR$Log2_64PP<-log2(UV_dmg_repair_mutR$Mean_64PP_FC)
UV_dmg_repair_mutR$Log2_CPD<-log2(UV_dmg_repair_mutR$Mean_CPD_FC)
mcols(UV_dmg_repair_mutR)


## Add UV ratio scores to master dataframe
UV_dmg_repair_mutR$UV64_ratio<-0
UV_dmg_repair_mutR$Log2_UV64_ratio<-0

UV_dmg_repair_mutR$UVcpd_ratio<-0
UV_dmg_repair_mutR$Log2_UVcpd_ratio<-0


UV_dmg_repair_mutR
NROW(UV_dmg_repair_mutR)


for (i in 1:NROW(UV_dmg_repair_mutR)){
  ID<-UV_dmg_repair_mutR$name[i]
  
  if (nrow(IP64_ratio_MutR[which(IP64_ratio_MutR$name==ID),])==1) {
    UV_dmg_repair_mutR$UV64_ratio[i]<-IP64_ratio_MutR[which(IP64_ratio_MutR$name==ID),11]
    UV_dmg_repair_mutR$Log2_UV64_ratio[i]<-IP64_ratio_MutR[which(IP64_ratio_MutR$name==ID),12]
  } else {
    UV_dmg_repair_mutR$UV64_ratio[i]<-NA
    UV_dmg_repair_mutR$Log2_UV64_ratio[i]<-NA
  }

  if (nrow(CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$name==ID),])==1) {
    UV_dmg_repair_mutR$UVcpd_ratio[i]<-CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$name==ID),11]
    UV_dmg_repair_mutR$Log2_UVcpd_ratio[i]<-CPD_dmg_repair_mutR[which(CPD_dmg_repair_mutR$name==ID),12]
  } else {
    UV_dmg_repair_mutR$UVcpd_ratio[i]<-NA
    UV_dmg_repair_mutR$Log2_UVcpd_ratio[i]<-NA
  }
  
  if (i==(NROW(UV_dmg_repair_mutR)/4)*1){
    print("25% done")
  }
  if (i==(NROW(UV_dmg_repair_mutR)/4)*2){
    print("50% done")
  }
  if (i==(NROW(UV_dmg_repair_mutR)/4)*3){
    print("75% done")
  }
  
}

UV_dmg_repair_mutR[15,]

Enh_UVdmg_repair_mutR<-UV_dmg_repair_mutR
Enh_UVdmg_repair_mutR<-as.data.frame(Enh_UVdmg_repair_mutR)
head(Enh_UVdmg_repair_mutR)

save(Enh_UVdmg_repair_mutR,file = paste0(out, "Enhancers_w_repair_mutR_UVratio_info.RData"))


### END ###

library(GenomicRanges)
library(rtracklayer)
library(EnsDb.Hsapiens.v75)


## Description: Measure NER rates within the gene body of protein coding genes using only XRseq signals that uniquely map to the genome.


## Strategy:
# Step 1: Retrieve list of all human protein coding genes and their genomic coordinates. 
# Step 2: Prepare for loop to process and bin XRseq data.
# Step 3: Calculate repair signal for both strands of a gene regardless of which side the reference gene is on. Use binning function
# Step 4: Calculate the average NER signal at each repair time point using both replicate XRseq files.
# (Measure Repair signal on both strands of a gene regardless of which side the reference gene is on.)
# Step 5: Average both replicates together.


out<-"output_directory/"


######
## Step 1: Retrieve list of all human protein coding genes and their genomic coordinates. 
EnsDb<-EnsDb.Hsapiens.v75
EnsDb_genes<-genes(EnsDb,columns = c(listColumns(EnsDb, "gene"), "entrezid"))
seqlevelsStyle(EnsDb_genes)<-"UCSC"
EnsDb_genes<-keepSeqlevels(EnsDb_genes, value = standardChromosomes(EnsDb_genes),pruning.mode = "coarse")
genome(EnsDb_genes)<-"hg19"
seqinfo(EnsDb_genes)
EnsDb_genes$seq_coord_system<-NULL
EnsDb_genes$symbol<-NULL
EnsDb_genes

# Keep genes with "protein_coding" bio_type descriptions
Prot_coding<-EnsDb_genes[which(EnsDb_genes$gene_biotype == "protein_coding"),] # Remove rows with no protein_id
Prot_coding<-sortSeqlevels(Prot_coding)
Prot_coding<-sort(Prot_coding)


######
## Step 2: Prepare for loop to process and bin XRseq data.

## XR64 repair data

out_64<-out
XR64_path<-"path_to_uniquely_mappable_XRseq_reads/"

# File names
file_names<-list.files(XR64_path) # 20 files
PatternList<-c("5min","20min","1h","2h","4h")


for (i in 1:5) {
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "Plus",value = T)
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "Plus",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "Minus",value = T)
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "Minus",value = T)
  
  ## Plus strand
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(XR64_path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep1_plus_bw)
  Prot_coding<-sort(Prot_coding)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep1_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_plus_Prot<-as.data.frame(Rep1_plus_Prot)
  print("Binning Rep1_plus_Prot complete")
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(XR64_path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep2_plus_bw)
  Prot_coding<-sort(Prot_coding)
  
  Rep2_plus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep2_plus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_plus_Prot<-as.data.frame(Rep2_plus_Prot)
  print("Binning Rep2_plus_Prot complete")
  
  # Combine both replicates and average them together
  Pooled_plus_Prot<-cbind(Rep1_plus_Prot, Rep2_plus_Prot[,10])
  Pooled_plus_Prot$R1_R2<-0
  Pooled_plus_Prot$R1_R2<-rowMeans(Pooled_plus_Prot[,c(10:11)], dims = 1)
  colnames(Pooled_plus_Prot)[c(10:12)]<-c("XR64_plus_Rep1","XR64_plus_Rep2","XR64_plus_mean")
  
  Pooled_plus_Prot<-makeGRangesFromDataFrame(Pooled_plus_Prot, keep.extra.columns = T)
  Pooled_plus_Prot<-sortSeqlevels(Pooled_plus_Prot)
  Pooled_plus_Prot<-sort(Pooled_plus_Prot)
  
  save(Pooled_plus_Prot,file = paste0(out_64, "XR64PP_PLUS_",  PatternList[i], "_Prot_genes.RData"))
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(XR64_path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep1_minus_bw)
  Prot_coding<-sort(Prot_coding)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep1_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep1_minus_Prot<-as.data.frame(Rep1_minus_Prot)
  print("Binning Rep1_minus_Prot complete")
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(XR64_path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep2_minus_bw)
  Prot_coding<-sort(Prot_coding)
  
  Rep2_minus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep2_minus_bw, varname = "XR64_average", na.rm = TRUE)
  Rep2_minus_Prot<-as.data.frame(Rep2_minus_Prot)
  print("Binning Rep2_minus_Prot complete")
  
  # Combine both replicates and average them together
  Pooled_minus_Prot<-cbind(Rep1_minus_Prot, Rep2_minus_Prot[,10])
  Pooled_minus_Prot$R1_R2<-0
  Pooled_minus_Prot$R1_R2<-rowMeans(Pooled_minus_Prot[,c(10:11)], dims = 1)
  colnames(Pooled_minus_Prot)[c(10:12)]<-c("XR64_minus_Rep1","XR64_minus_Rep2","XR64_minus_mean")
  
  Pooled_minus_Prot<-makeGRangesFromDataFrame(Pooled_minus_Prot, keep.extra.columns = T)
  Pooled_minus_Prot<-sortSeqlevels(Pooled_minus_Prot)
  Pooled_minus_Prot<-sort(Pooled_minus_Prot)
  
  save(Pooled_minus_Prot,file = paste0(out_64, "XR64PP_MINUS_",  PatternList[i], "_Prot_genes.RData"))
  
}


## XR CPD repair data

out_CPD<-out
XRCPD_path<-"path_to_uniquely_mappable_XRseq_reads/"

# File names
file_names<-list.files(XRCPD_path)
PatternList<-c("1h","4h","8h","16h","1d","2d")


for (i in 1:6) {
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "Plus",value = T)
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "Plus",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "Minus",value = T)
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "Minus",value = T)
  
  ## Plus strand
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep1_plus),as='RleList')
  Rep1_plus_bw[Rep1_plus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep1_plus_bw)
  Prot_coding<-sort(Prot_coding)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_plus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep1_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_plus_Prot<-as.data.frame(Rep1_plus_Prot)
  print("Binning Rep1_plus_Prot complete")
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep2_plus),as='RleList')
  Rep2_plus_bw[Rep2_plus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep2_plus_bw)
  Prot_coding<-sort(Prot_coding)
  
  Rep2_plus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep2_plus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_plus_Prot<-as.data.frame(Rep2_plus_Prot)
  print("Binning Rep2_plus_Prot complete")
  
  # Combine both replicates and average them together
  Pooled_plus_Prot<-cbind(Rep1_plus_Prot, Rep2_plus_Prot[,10])
  Pooled_plus_Prot$R1_R2<-0
  Pooled_plus_Prot$R1_R2<-rowMeans(Pooled_plus_Prot[,c(10:11)], dims = 1)
  colnames(Pooled_plus_Prot)[c(10:12)]<-c("XRCPD_plus_Rep1","XRCPD_plus_Rep2","XRCPD_plus_mean")
  
  Pooled_plus_Prot<-makeGRangesFromDataFrame(Pooled_plus_Prot, keep.extra.columns = T)
  Pooled_plus_Prot<-sortSeqlevels(Pooled_plus_Prot)
  Pooled_plus_Prot<-sort(Pooled_plus_Prot)
  
  save(Pooled_plus_Prot,file = paste0(out_CPD, "XRCPD_PLUS_",  PatternList[i], "_Prot_genes.RData"))
  
  ## Minus strand
  # Rep1
  Rep1_minus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep1_minus),as='RleList')
  Rep1_minus_bw[Rep1_minus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep1_minus_bw)
  Prot_coding<-sort(Prot_coding)
  
  # Bin XRseq signal to genes on plus strand
  Rep1_minus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep1_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep1_minus_Prot<-as.data.frame(Rep1_minus_Prot)
  print("Binning Rep1_minus_Prot complete")
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(XRCPD_path,"/",Rep2_minus),as='RleList')
  Rep2_minus_bw[Rep2_minus_bw == 0]<-NA
  
  seqlevels(Prot_coding)<-seqlevels(Rep2_minus_bw)
  Prot_coding<-sort(Prot_coding)
  
  Rep2_minus_Prot<-binnedAverage(bins = Prot_coding, numvar = Rep2_minus_bw, varname = "XRCPD_average", na.rm = TRUE)
  Rep2_minus_Prot<-as.data.frame(Rep2_minus_Prot)
  print("Binning Rep2_minus_Prot complete")
  
  # Combine both replicates and average them together
  Pooled_minus_Prot<-cbind(Rep1_minus_Prot, Rep2_minus_Prot[,10])
  Pooled_minus_Prot$R1_R2<-0
  Pooled_minus_Prot$R1_R2<-rowMeans(Pooled_minus_Prot[,c(10:11)], dims = 1)
  colnames(Pooled_minus_Prot)[c(10:12)]<-c("XRCPD_minus_Rep1","XRCPD_minus_Rep2","XRCPD_minus_mean")
  
  Pooled_minus_Prot<-makeGRangesFromDataFrame(Pooled_minus_Prot, keep.extra.columns = T)
  Pooled_minus_Prot<-sortSeqlevels(Pooled_minus_Prot)
  Pooled_minus_Prot<-sort(Pooled_minus_Prot)
  
  save(Pooled_minus_Prot,file = paste0(out_CPD, "XRCPD_MINUS_",  PatternList[i], "_Prot_genes.RData"))
  
}


## END ##
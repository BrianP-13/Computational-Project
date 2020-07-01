library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


## Description: This script measures the C>T (and G>A) mutation rates across the genome, binned at 100kb and 1Mb size.


out<-"output_directory/"



## Import cutaneous SSM dataset
load(file = "file_path/MELA_AU_cutaneous_SSM.RData")
SSM<-SSM[,c(3,4,5,7,9)]
SSM<-makeGRangesFromDataFrame(SSM, keep.extra.columns = T,
                              start.field = "chromosome_start",
                              end.field = "chromosome_end")
seqlevelsStyle(SSM)<-"UCSC"
SSM<-sortSeqlevels(SSM)
SSM<-sort(SSM)


######
## Measuring mutation rates in 1Mb genome

# create 1Mb binned genome
gen <- BSgenome.Hsapiens.UCSC.hg19
si.gen <- seqinfo(gen)
si<-keepSeqlevels(si.gen, value = standardChromosomes(si.gen),pruning.mode = "coarse")
si<-seqlengths(si)
SBS_1Mb<-tileGenome(si,tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)


# Create columns with mutation counts to add to GRanges
ColNames<-c('G_to_A_count','G_to_T_count','G_to_C_count','C_to_T_count','C_to_A_count','C_to_G_count','T_to_A_count','T_to_G_count','T_to_C_count','A_to_T_count','A_to_C_count','A_to_G_count')
df<-as.data.frame(matrix(ncol =12 ,nrow =3114 ))
colnames(df)<-ColNames

mcols(SBS_1Mb)<-cbind(mcols(SBS_1Mb), df)

# Create list of reference nt> mutated nt base pairs
Mut_types<-list(c("G","A"),c("G","T"),c("G","C"),
                c("C","T"),c("C","A"),c("C","G"),
                c("T","A"),c("T","G"),c("T","C"),
                c("A","T"),c("A","C"),c("A","G"))


## Measure mutated base counts per genomic region
for (j in 1:12) { 
  mutN<-unlist(Mut_types[j])
  N_to_mutN<-subset(SSM,reference_genome_allele==mutN[1])
  N_to_mutN<-subset(N_to_mutN, mutated_to_allele==mutN[2])
  
  for (i in 1:3114) { 
    hits<-findOverlaps(query = N_to_mutN,subject = SBS_1Mb[i])
    mcols(SBS_1Mb)[i,j]<-length(hits)
  }
}

## Next normalize counts by reference  nucleotide  content
ColNames2<-c('Normalized_G_to_A_count','Normalized_G_to_T_count','Normalized_G_to_C_count','Normalized_C_to_T_count','Normalized_C_to_A_count','Normalized_C_to_G_count','Normalized_T_to_A_count','Normalized_T_to_G_count','Normalized_T_to_C_count','Normalized_A_to_T_count','Normalized_A_to_C_count','Normalized_A_to_G_count')
df2<-as.data.frame(mcols(SBS_1Mb))
colnames(df2)<-ColNames2

# Get DNA sequence info for human genome at 1Mb DNA lengths
BSgenome<-BSgenome.Hsapiens.UCSC.hg19
DNA<-getSeq(x = BSgenome,names=SBS_1Mb)

# Calculate individual nt content per genomic region
G_content<-letterFrequency(x = DNA, letters = "G")
C_content<-letterFrequency(x = DNA, letters = "C")
T_content<-letterFrequency(x = DNA, letters = "T")
A_content<-letterFrequency(x = DNA, letters = "A")

# Normalize mutated base counts by reference base content per region
df2[,c(1:3)]<-df2[,c(1:3)]/G_content
df2[,c(4:6)]<-df2[,c(4:6)]/C_content
df2[,c(7:9)]<-df2[,c(7:9)]/T_content
df2[,c(10:12)]<-df2[,c(10:12)]/A_content


# Combine normalized mutation rates to master table
mcols(SBS_1Mb)<-cbind(mcols(SBS_1Mb), df2)

## Save work

save(SBS_1Mb, file = paste0(out,'All_SSM_rates_1Mb_hg19.RData'))
#write.csv(x = as.data.frame(SBS_1Mb), file = paste0(path,'Mutation rate files/', 'All_SSM_rates_1Mb_hg19.csv'))


######
## Measuring mutation rates in 100kb genome

# create 100kb binned genome
gen <- BSgenome.Hsapiens.UCSC.hg19
si.gen <- seqinfo(gen)
si<-keepSeqlevels(si.gen, value = standardChromosomes(si.gen),pruning.mode = "coarse")
si<-seqlengths(si)
SBS_100kb<-tileGenome(si,tilewidth = 1e5, cut.last.tile.in.chrom = TRUE)

# Create columns with mutation counts to add to GRanges
ColNames<-c('G_to_A_count','G_to_T_count','G_to_C_count','C_to_T_count','C_to_A_count','C_to_G_count','T_to_A_count','T_to_G_count','T_to_C_count','A_to_T_count','A_to_C_count','A_to_G_count')
df<-as.data.frame(matrix(ncol =12 ,nrow =30971 ))
colnames(df)<-ColNames

mcols(SBS_100kb)<-cbind(mcols(SBS_100kb), df)

# Create list of reference nt> mutated nt base pairs
Mut_types<-list(c("G","A"),c("G","T"),c("G","C"),
                c("C","T"),c("C","A"),c("C","G"),
                c("T","A"),c("T","G"),c("T","C"),
                c("A","T"),c("A","C"),c("A","G"))

## Measure mutated base counts per genomic region
for (j in 1:12) { 
  mutN<-unlist(Mut_types[j])
  N_to_mutN<-subset(SSM,reference_genome_allele==mutN[1])
  N_to_mutN<-subset(N_to_mutN, mutated_to_allele==mutN[2])
  
  for (i in 1:30971) { 
    hits<-findOverlaps(query = N_to_mutN,subject = SBS_100kb[i])
    mcols(SBS_100kb)[i,j]<-length(hits)
  }
  print(paste0("Done with ",mutN[1],">",mutN[2]))
}

## Next normalize counts by reference  nucleotide  content
ColNames2<-c('Normalized_G_to_A_count','Normalized_G_to_T_count','Normalized_G_to_C_count','Normalized_C_to_T_count','Normalized_C_to_A_count','Normalized_C_to_G_count','Normalized_T_to_A_count','Normalized_T_to_G_count','Normalized_T_to_C_count','Normalized_A_to_T_count','Normalized_A_to_C_count','Normalized_A_to_G_count')
df2<-as.data.frame(mcols(SBS_100kb))
colnames(df2)<-ColNames2

# Get DNA sequence info for human genome at 1Mb DNA lengths
BSgenome<-BSgenome.Hsapiens.UCSC.hg19
DNA<-getSeq(x = BSgenome,names=SBS_100kb)

# Calculate individual nt content per genomic region
G_content<-letterFrequency(x = DNA, letters = "G")
C_content<-letterFrequency(x = DNA, letters = "C")
T_content<-letterFrequency(x = DNA, letters = "T")
A_content<-letterFrequency(x = DNA, letters = "A")

# Normalize mutated base counts by reference base content per region
df2[,c(1:3)]<-df2[,c(1:3)]/G_content
df2[,c(4:6)]<-df2[,c(4:6)]/C_content
df2[,c(7:9)]<-df2[,c(7:9)]/T_content
df2[,c(10:12)]<-df2[,c(10:12)]/A_content

# Combine normalized mutation rates to master table
mcols(SBS_100kb)<-cbind(mcols(SBS_100kb), df2)

## Save work
save(SBS_100kb, file = "file_path/All_SSM_rates_100kb_hg19.RData")



## END ##
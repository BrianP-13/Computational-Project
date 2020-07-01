library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)


### Processing script to subset single base mutations from cutaneous melanoma cancer samples (duplicate rows removed) ###


######
## Step 1: Subset SSM to include only cutaneous melanoma samples

# Import SSM file
path<-'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/Genomics Data/Cancer mutation data/Hayward et al 2017 MELA-WGS/Mutation data files/'
SSM<-read.csv(paste0(path, 'SSM_essential_info.MELA_AU.csv.gz')) # actual SSM from study without non essential columns to reduce size of file
SSM<-SSM[,c(2:18)]
head(SSM)

## Need to subset table to only include samples from cutaneous melanoma cancers
path2<-'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/Genomics Data/Cancer mutation data/Hayward et al 2017 MELA-WGS/Donor ID files/'
Donor_ID<-read.csv(paste0(path2, "Donor_ID_list_cutaneous_MELA_with_UV_sig.csv"), header = T)
Donor_ID<-Donor_ID[,c(2:3)]
Donor_ID[,1]<-levels(droplevels(Donor_ID[,1]))
Donor_ID[,2]<-levels(droplevels(Donor_ID[,2]))
head(Donor_ID)

## Get indeces for rows in SSM that match the donor_ID list
Rows<-list()
for (i in 1:136) { 
  indices<-which(SSM$icgc_donor_id == Donor_ID[i,1])
  Rows<-append(Rows, indices)
}
Rows<-unlist(Rows)
SSM<-SSM[Rows,]

#SSM<-read.csv(paste0(path, 'simple_somatic_mutations.cutaneous.MELA_AU.csv.gz')) # contains duplicates

## Keep only relevant information with mutation data
head(SSM)

SSM<-SSM[,c(2,1,9,10,11,14,15,16,17)] # rows to keep
head(SSM)

# De-duplicate rows
SSM<-distinct(SSM)

save(SSM, file = paste0(path,'MELA_AU_cutaneous_SSM.RData'),compress = "gzip")


# to load mutation file:

load(file = 'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/Genomics Data/Cancer mutation data/Hayward et al 2017 MELA-WGS/Mutation data files/MELA_AU_cutaneous_SSM.RData')
SSM


#length(unique(SSM$icgc_mutation_id)) # number of unique mut IDs = 21304761 (SAME!), number of rows = 22961363
#path<-'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/Genomics Data/Cancer Mutations/MELA-AU/'
#write.csv( SSM, paste0(path, "MELA_AU_unique_SSM_essential_info.csv"))




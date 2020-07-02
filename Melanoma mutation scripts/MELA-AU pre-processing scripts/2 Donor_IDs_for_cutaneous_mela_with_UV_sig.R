

# Hayward et al, 2017 , DOI: 10.1038/nature22071
path<-'D:/Brian/Google Drive/Important Documents/Stanford University/LAB/Computational Project/Genomics Data/Cancer mutation data/Hayward et al 2017 MELA-WGS/Donor ID files/'

# Donor file with melanoma subtypes
Subtypes<-read.csv(paste0(path, 'Donor_melanoma_subtypes.csv'), header = T)
head(Subtypes)

# Subset Only rows with cutaneous melanoma samples and have UV signature, should have 136 tumor samples
Cutaneous<-Subtypes[which(Subtypes$subtype == "Cutaneous"),]
Cutaneous_sig7<-Cutaneous[which(Cutaneous$Major_signature..UV.or.non.UV. =="UV"),]
head(Cutaneous_sig7)

Donor_List<-unlist( as.list(Cutaneous_sig7[,1]))
Donor_List<-as.character(Donor_List)

donorL<-list()
for (i in 1:length(Donor_List)) { 
  names<-unlist(strsplit(Donor_List[i], "_"))
  names<-paste0(names[1],'-',names[2])
  donorL<-append(donorL, names)
}
donorL<-unlist(donorL)

# Donor file
Donor<-read.table(file = paste0(path, 'donor.MELA-AU.tsv.gz'),header = T,sep = "\t")
Donor<-Donor[,c(1,4)]
Donor[,1]<-levels(droplevels(Donor[,1]))
Donor[,2]<-levels(droplevels(Donor[,2]))

# Next, match the names from "Subtype" to the donor ids in the Donor df

IDs<-list()
for (i in 1:length(Donor_List)) { 
  name<-Donor[which(Donor$submitted_donor_id== donorL[i]),]
  name$icgc_donor_id
  IDs<-append(IDs,name$icgc_donor_id)
}
IDs<-unlist(IDs)

Names<-as.data.frame(matrix(ncol = 2,nrow = 136))
colnames(Names)<-c("icgc_donor_id","submitted_donor_id")
Names[,1]<-IDs
Names[,2]<-donorL

write.csv(Names, paste0(path, "Donor_ID_list_cutaneous_MELA_with_UV_sig.csv"))


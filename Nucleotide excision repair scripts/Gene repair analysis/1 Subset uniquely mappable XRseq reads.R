library(GenomicRanges)
library(rtracklayer)


## Description: This script subsets only XRseq signals that uniquely map to the genome within 24 bp bins


## Strategy:
# Step 1: Import the ENCODE_24mer_mappability file to get coordinates with unique 24mer sites in the genome
# Step 2: Set up forloop to import each XRseq replicate file for both strands
# Step 3: Save new XRseq with unique regions


out<-"output_directory/"


######
## Step 1: Import the ENCODE_24mer_mappability file to get coordinates with unique 24mer sites in the genome

# Import ENCODE24mer_mappability_score1 file with coordinates for unique regions
load(file = "file_path/ENCODE24mer_mappability_score1.RData")
ENCODE_24_mappability


######
## 64PP repair


# paths to 64PP XRseq bw files
path<-file.path("path_to_XRseq_files/")

# File names
file_names<-list.files("path_to_XRseq_files/") # 20 files
PatternList<-c("5min","20min","1h","2h","4h")


for (i in 1:5) {
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "PLUS",value = T)
  
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "PLUS",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "MINUS",value = T)
  
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "MINUS",value = T)
  
  # Plus strand Rep1 and Rep2
  
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus)) # import as GRanges object class
  Rep1_plus_bw<-subsetByOverlaps(x = Rep1_plus_bw,ranges = ENCODE_24_mappability,minoverlap = 1,type = "start")
  export.bw(object = Rep1_plus_bw, con = paste0(out,'NHF164_',PatternList[i],'_Rep1_Plus_strand_unique.bw'))
  print(paste0("Done with XR64 ", PatternList[i],' Rep1 PLUS'))
  rm(Rep1_plus_bw)
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus))
  Rep2_plus_bw<-subsetByOverlaps(x = Rep2_plus_bw,ranges = ENCODE_24_mappability, minoverlap = 1,type = "start")
  export.bw(object = Rep2_plus_bw, con = paste0(out,'NHF164_',PatternList[i],'_Rep2_Plus_strand_unique.bw'))
  print(paste0("Done with XR64 ", PatternList[i],' Rep2 PLUS'))
  rm(Rep2_plus_bw)
  
  # Minus strand Rep1 and Rep2
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus)) # import as GRanges object class
  Rep1_minus_bw<-subsetByOverlaps(x = Rep1_minus_bw,ranges = ENCODE_24_mappability,minoverlap = 1,type = "start")
  export.bw(object = Rep1_minus_bw, con = paste0(out,'NHF164_',PatternList[i],'_Rep1_Minus_strand_unique.bw'))
  print(paste0("Done with XR64 ", PatternList[i],' Rep1 MINUS'))
  rm(Rep1_minus_bw)
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus))
  Rep2_minus_bw<-subsetByOverlaps(x = Rep2_minus_bw,ranges = ENCODE_24_mappability, minoverlap = 1,type = "start")
  export.bw(object = Rep2_minus_bw, con = paste0(out,'NHF164_',PatternList[i],'_Rep2_Minus_strand_unique.bw'))
  print(paste0("Done with XR64 ", PatternList[i],' Rep2 MINUS'))
  rm(Rep2_minus_bw)
  
}


######
## CPD repair

"path_to_XRseq_files/"
# paths to CPD XRseq bw files
path<-file.path("path_to_XRseq_files/")
file_names<-list.files("path_to_XRseq_files/")
PatternList<-c("1h","4h","8h","16h","1d","2d")

for (i in 1:6) {
  Rep1_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_plus<-grep(x = Rep1_plus,pattern = "Rep1",value = T)
  Rep1_plus<-grep(x = Rep1_plus, pattern = "PLUS",value = T)
  
  Rep2_plus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_plus<-grep(x = Rep2_plus,pattern = "Rep2",value = T)
  Rep2_plus<-grep(x = Rep2_plus, pattern = "PLUS",value = T)
  
  Rep1_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep1_minus<-grep(x = Rep1_minus,pattern = "Rep1",value = T)
  Rep1_minus<-grep(x = Rep1_minus, pattern = "MINUS",value = T)
  
  Rep2_minus<-grep(x = file_names,pattern = PatternList[i],value = T)
  Rep2_minus<-grep(x = Rep2_minus,pattern = "Rep2",value = T)
  Rep2_minus<-grep(x = Rep2_minus, pattern = "MINUS",value = T)
  
  # Plus strand Rep1 and Rep2
  
  # Rep1
  Rep1_plus_bw<-import.bw(con = paste0(path,"/",Rep1_plus)) # import as GRanges object class
  Rep1_plus_bw<-subsetByOverlaps(x = Rep1_plus_bw,ranges = ENCODE_24_mappability,minoverlap = 1,type = "start")
  export.bw(object = Rep1_plus_bw, con = paste0(out,'NHF1CPD_',PatternList[i],'_Rep1_Plus_strand_unique.bw'))
  print(paste0("Done with XR CPD ", PatternList[i],' Rep1 PLUS'))
  rm(Rep1_plus_bw)
  
  # Rep2
  Rep2_plus_bw<-import.bw(con = paste0(path,"/",Rep2_plus))
  Rep2_plus_bw<-subsetByOverlaps(x = Rep2_plus_bw,ranges = ENCODE_24_mappability, minoverlap = 1,type = "start")
  export.bw(object = Rep2_plus_bw, con = paste0(out,'NHF1CPD_',PatternList[i],'_Rep2_Plus_strand_unique.bw'))
  print(paste0("Done with XR CPD ", PatternList[i],' Rep2 PLUS'))
  rm(Rep2_plus_bw)
  
  # Minus strand Rep1 and Rep2
  Rep1_minus_bw<-import.bw(con = paste0(path,"/",Rep1_minus)) # import as GRanges object class
  Rep1_minus_bw<-subsetByOverlaps(x = Rep1_minus_bw,ranges = ENCODE_24_mappability,minoverlap = 1,type = "start")
  export.bw(object = Rep1_minus_bw, con = paste0(out,'NHF1CPD_',PatternList[i],'_Rep1_Minus_strand_unique.bw'))
  print(paste0("Done with XR CPD ", PatternList[i],' Rep1 MINUS'))
  rm(Rep1_minus_bw)
  
  # Rep2
  Rep2_minus_bw<-import.bw(con = paste0(path,"/",Rep2_minus))
  Rep2_minus_bw<-subsetByOverlaps(x = Rep2_minus_bw,ranges = ENCODE_24_mappability, minoverlap = 1,type = "start")
  export.bw(object = Rep2_minus_bw, con = paste0(out,'NHF1CPD_',PatternList[i],'_Rep2_Minus_strand_unique.bw'))
  print(paste0("Done with XR CPD ", PatternList[i],' Rep2 MINUS'))
  rm(Rep2_minus_bw)
  
}


## END ##

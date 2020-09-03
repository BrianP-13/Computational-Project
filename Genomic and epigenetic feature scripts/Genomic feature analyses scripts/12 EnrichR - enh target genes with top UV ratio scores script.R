library(enrichR)


## Description: Plot with all repeat class enrichment results from all regions of interest.


out<-"output_directory/"


dbs <- c("GO_Biological_Process_2018","GO_Molecular_Function_2018","GO_Cellular_Component_2018","KEGG_2019_Human")


## Enh-Genes in top/bottom 10% 6-4PP ratio regions
T10_64_genes<-read.csv(file = "file_path/Enh_gene_list_in_Top10_64PP_ratio.csv")
B10_64_genes<-read.csv(file = "file_path/Enh_gene_list_in_Bot10_64PP_ratio.csv")

## Enh-Genes in top/bottom 10% CPD ratio regions
T10_CPD_genes<-read.csv(file = "file_path/Enh_gene_list_in_Top10_CPD_ratio.csv")
B10_CPD_genes<-read.csv(file = "file_path/Enh_gene_list_in_Bot10_CPD_ratio.csv")


# List of Enhancer-target genes in top and bottom regions, file name components
Gene_list_T10<-list(T10_64_genes,T10_CPD_genes)
Gene_list_B10<-list(B10_64_genes,B10_CPD_genes)

Path_out<-c("IP64/","CPD/")
File_names_T10<-c("T10_64_ratio_enh_genes","T10_CPD_ratio_enh_genes")
File_names_B10<-c("B10_64_ratio_enh_genes","B10_CPD_ratio_enh_genes")

for (i in 1:2) {
  ## Genes in top regions
  List_T10<-Gene_list_T10[[i]]
  List_T10<-as.character(List_T10$hgnc_symbol)
  List_T10<-na.omit(List_T10)
  List_T10<-List_T10[which(List_T10 != "")]
  
  Gene_Enrch_T10<-enrichr(genes = List_T10,databases =  dbs)
  
  # GO Biological Process
  Gene_Enrch_T10_GO_BP<-Gene_Enrch_T10$GO_Biological_Process_2018
  Gene_Enrch_T10_GO_BP$Overlap<-paste("'",Gene_Enrch_T10_GO_BP$Overlap)
  write.csv(Gene_Enrch_T10_GO_BP, file = paste0(out,Path_out[i],File_names_T10[i],"_GO_BP.csv"))
  
  # GO Molecular Function
  Gene_Enrch_T10_GO_MF<-Gene_Enrch_T10$GO_Molecular_Function_2018
  Gene_Enrch_T10_GO_MF$Overlap<-paste("'",Gene_Enrch_T10_GO_MF$Overlap)
  write.csv(Gene_Enrch_T10_GO_MF, file = paste0(out,Path_out[i],File_names_T10[i],"_GO_MF.csv"))
  # GO Cellular Component
  Gene_Enrch_T10_GO_CC<-Gene_Enrch_T10$GO_Cellular_Component_2018
  Gene_Enrch_T10_GO_CC$Overlap<-paste("'",Gene_Enrch_T10_GO_CC$Overlap)
  write.csv(Gene_Enrch_T10_GO_CC, file = paste0(out,Path_out[i],File_names_T10[i],"_GO_CC.csv"))
  
  # KEGG
  Gene_Enrch_T10_KEGG<-Gene_Enrch_T10$KEGG_2019_Human
  Gene_Enrch_T10_KEGG$Overlap<-paste("'",Gene_Enrch_T10_KEGG$Overlap)
  write.csv(Gene_Enrch_T10_KEGG, file = paste0(out,Path_out[i],File_names_T10[i],"_GO_KEGG.csv"))
  
  print(paste0("Done with ",File_names_T10[i]))
  
  ## Genes in bottom regions
  List_B10<-Gene_list_B10[[i]]
  List_B10<-as.character(List_B10$hgnc_symbol)
  List_B10<-na.omit(List_B10)
  List_B10<-List_B10[which(List_B10 != "")]
  
  Gene_Enrch_B10<-enrichr(genes = List_B10,databases =  dbs)
  
  # GO Biological Process
  Gene_Enrch_B10_GO_BP<-Gene_Enrch_B10$GO_Biological_Process_2018
  Gene_Enrch_B10_GO_BP$Overlap<-paste("'",Gene_Enrch_B10_GO_BP$Overlap)
  write.csv(Gene_Enrch_B10_GO_BP, file = paste0(out,Path_out[i],File_names_B10[i],"_GO_BP.csv"))
  
  # GO Molecular Function
  Gene_Enrch_B10_GO_MF<-Gene_Enrch_B10$GO_Molecular_Function_2018
  Gene_Enrch_B10_GO_MF$Overlap<-paste("'",Gene_Enrch_B10_GO_MF$Overlap)
  write.csv(Gene_Enrch_B10_GO_MF, file = paste0(out,Path_out[i],File_names_B10[i],"_GO_MF.csv"))
  
  # GO Cellular Component
  Gene_Enrch_B10_GO_CC<-Gene_Enrch_B10$GO_Cellular_Component_2018
  Gene_Enrch_B10_GO_CC$Overlap<-paste("'",Gene_Enrch_B10_GO_CC$Overlap)
  write.csv(Gene_Enrch_B10_GO_CC, file = paste0(out,Path_out[i],File_names_B10[i],"_GO_CC.csv"))
  
  # KEGG
  Gene_Enrch_B10_KEGG<-Gene_Enrch_B10$KEGG_2019_Human
  Gene_Enrch_B10_KEGG$Overlap<-paste("'",Gene_Enrch_B10_KEGG$Overlap)
  write.csv(Gene_Enrch_B10_KEGG, file = paste0(out,Path_out[i],File_names_B10[i],"_GO_KEGG.csv"))
  
  print(paste0("Done with ",File_names_B10[i]))
  
}



## END ##
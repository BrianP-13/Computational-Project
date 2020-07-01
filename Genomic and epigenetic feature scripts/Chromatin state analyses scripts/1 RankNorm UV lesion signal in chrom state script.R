library(GenomicRanges)
library(ggplot2)

# Description: This script takes the UV lesion fold change signal binned to IMR90 15 chromatin states and rank normalizes them to measure the relative differences in susceptibility


out<-"output_directory/"

# rank-based inverse normalization normal transformation function
rankitNormalize_vector <- function(x) {
  
  stopifnot(is.numeric(x))
  x <- qnorm((rank(x) - 0.5) / length(x))
  return(x)
}

rankitNormalize <- function(x, IND = 1) {
  
  # Normalizes rows (IND = 1) or columns (IND = 2)
  # to a quantile standard normalization (i.e. rankit)
  
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(x))
  
  rowNames <- rownames(x)
  colNames <- colnames(x)
  
  x <- apply(x, IND, rankitNormalize_vector)
  if(IND == 1)
    x <- t(x)
  
  rownames(x) <- rowNames
  colnames(x) <- colNames
  
  return(x)
  
}


## Combine binned UV lesion signals together in one data frame

# Load completed CPD Datatable here:
load(file ="file_path/CPD_binned_to_chrom15.RData")
dfCPD15<-as.data.frame(CPD_15)
colnames(dfCPD15)[6]<-"State"
colnames(dfCPD15)[7]<-"CPD_mean"

# Load completed 6-4PP Datatable here:
load(file ="file_path/IP64_binned_to_chrom15.RData")
df64<-as.data.frame(IP64_15)

Lesions_in_ChromState<-cbind(dfCPD15, df64$IP64_mean)
colnames(Lesions_in_ChromState)[8]<-"IP64_mean"
Lesions_in_ChromState<-na.omit(Lesions_in_ChromState)

# Apply rank normalization function
CPD_64PP<-Lesions_in_ChromState[,c(7,8)]
CPD_64PP<-as.matrix(CPD_64PP)

ranked_UV_lesions<-rankitNormalize(x= CPD_64PP, IND = 2)
ranked_UV_lesions<-as.data.frame(ranked_UV_lesions)


# Add rank normalized values back to main dataset
Lesions_in_ChromState<-cbind(Lesions_in_ChromState, ranked_UV_lesions)
colnames(Lesions_in_ChromState)[c(9,10)]<-c("ranked_CPD","ranked_64PP")

hist(Lesions_in_ChromState$CPD_mean, breaks = 100,
     main = "CPD 1Mb bin count vs signal",
     xlab = "CPD FC signal",
     ylab = "Bin counts")

hist(Lesions_in_ChromState$ranked_CPD, breaks = 100,
     main = "Ranked CPD bin count vs signal",
     xlab = "Ranked CPD signal",
     ylab = "Bin counts")

hist(Lesions_in_ChromState$IP64_mean, breaks = 100,
     main = "64PP 1Mb bin count vs signal",
     xlab = "64PP FC signal",
     ylab = "Bin counts")

hist(Lesions_in_ChromState$ranked_64PP, breaks = 100,
     main = "Ranked 64PP bin count vs signal",
     xlab = "Ranked UV signal",
     ylab = "Bin counts")

hist(Lesions_in_ChromState$ranked_64PP, breaks = 100,
     main = "Rank inverse normal transformed signal",
     xlab = "Ranked UV signal",
     ylab = "Bin counts")


# calculate distributions of rank differences in each chromatin state
Ranked_diff<-Lesions_in_ChromState[,c(1:6,9,10)]
Ranked_diff$Difference<-0
Ranked_diff$Difference<-Lesions$ranked_64PP-Lesions$ranked_CPD
hist(Ranked_diff$Difference)

save(Ranked_diff, file = paste0(out, "Ranked_UV_lesions_chrom_states.RData"))


## Ranked signal distribution in chromatin states

# 1. Ranked CPD signal

Box<-boxplot.stats(Ranked_diff$ranked_CPD)

# Measure data distribution using quantile function
CPDdf = as.data.frame(matrix(nrow=16, ncol=7))
colnames(CPDdf) = c("states","min","q25","median","q75","max","mean")

CPDdf$states[16] = "Whole genome"
CPDdf$min[16] = Box$stats[1] # min
CPDdf$q25[16] = Box$stats[2] # q25
CPDdf$median[16] = Box$stats[3] # median
CPDdf$q75[16] = Box$stats[4] # q75
CPDdf$max[16] = Box$stats[5] # max
CPDdf$mean[16] = mean(Ranked_diff$ranked_CPD)

for (i in 1:15) {
  state<-Ranked_diff
  state<-state[which(state$State == chrom_state[i]),]
  
  box_stats<-boxplot.stats(state$ranked_CPD)
  
  CPDdf$states[i] = chrom_state[i]
  CPDdf$min[i] = box_stats$stats[1]
  CPDdf$q25[i] = box_stats$stats[2]
  CPDdf$median[i] = box_stats$stats[3]
  CPDdf$q75[i] = box_stats$stats[4]
  CPDdf$max[i] = box_stats$stats[5]
  CPDdf$mean[i] = mean(state$ranked_CPD)
}

save(CPDdf, file = paste0(out, "Ranked_CPD_boxplot_stats_for_chrom_states.RData"))

# 2. Ranked 64PP signal

Box<-boxplot.stats(Ranked_diff$ranked_64PP)

# Measure data distribution using quantile function
IP64df = as.data.frame(matrix(nrow=16, ncol=7))
colnames(IP64df) = c("states","min","q25","median","q75","max", "mean")

IP64df$states[16] = "Whole genome"
IP64df$min[16] = Box$stats[1] # min
IP64df$q25[16] = Box$stats[2] # q25
IP64df$median[16] = Box$stats[3] # median
IP64df$q75[16] = Box$stats[4] # q75
IP64df$max[16] = Box$stats[5] # max
IP64df$mean[16] = mean(Ranked_diff$ranked_64PP)

for (i in 1:15) {
  state<-Ranked_diff
  state<-state[which(state$State == chrom_state[i]),]
  
  box_stats<-boxplot.stats(state$ranked_64PP)
  
  IP64df$states[i] = chrom_state[i]
  IP64df$min[i] = box_stats$stats[1]
  IP64df$q25[i] = box_stats$stats[2]
  IP64df$median[i] = box_stats$stats[3]
  IP64df$q75[i] = box_stats$stats[4]
  IP64df$max[i] = box_stats$stats[5]
  IP64df$mean[i] = mean(state$ranked_64PP)
}

save(IP64df, file = paste0(out, "Ranked_64PP_boxplot_stats_for_chrom_states.RData"))

# 3. Difference

Box<-boxplot.stats(Ranked_diff$Difference)

Diffdf = as.data.frame(matrix(nrow=16, ncol=7))
colnames(Diffdf) = c("states","min","q25","median","q75","max","mean")

Diffdf$states[16] = "Whole genome"
Diffdf$min[16] = Box$stats[1] # min
Diffdf$q25[16] = Box$stats[2] # q25
Diffdf$median[16] = Box$stats[3] # median
Diffdf$q75[16] = Box$stats[4] # q75
Diffdf$max[16] = Box$stats[5] # max
Diffdf$mean[16] = mean(Ranked_diff$Difference)
Diffdf

for (i in 1:15) {
  state<-Ranked_diff
  state<-state[which(state$State == chrom_state[i]),]
  
  box_stats<-boxplot.stats(state$Difference)
  
  Diffdf$states[i] = chrom_state[i]
  Diffdf$min[i] = box_stats$stats[1]
  Diffdf$q25[i] = box_stats$stats[2]
  Diffdf$median[i] = box_stats$stats[3]
  Diffdf$q75[i] = box_stats$stats[4]
  Diffdf$max[i] = box_stats$stats[5]
  Diffdf$mean[i] = mean(state$Difference)
}

save(Diffdf, file = paste0(out, "Ranked_64PP_CPD_difference_boxplot_stats_for_chrom_states.RData"))


## END ##

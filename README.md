# UV Lesion Script organization
Date: September 3, 2020

## UV Lesion Script organization

Binning UV lesion signal scripts
1.	UV lesion binning (100kb and 1Mb) script
2.	Rank normalization of UV lesion signal script
3.	UV lesion binning per chromatin state script
4.	UV lesion binning per gene
5.	UV lesion binning per enhancer

UV lesion sequencing analyses
1.	IP 6-4PP replicate correlation script
2.	Di-pyrimidine frequency in raw sequencing reads script
3.	Di-pyrimidine frequency boxplot script

## Genomic and Epigenetic Feature Script organization
Genomic feature analyses
1.	RepeatMasker pre-processing script
2.	Enhancer target gene for IMR90 dataset EnhancerAtlas
3.	Genic feature enrichment in top UV susc regions scripts
4.	Enhancer element enrichment in top UV susc regions script
5.	Super enhancer enrichment in top UV susc regions script
6.	Genic feature and enhancer element enrichment boxplot script
7.	Repeat class enrichment in top UV susc regions script
8.	Repeat class boxplot script
9.	EnrichR - genes in top regions script
10.	EnrichR - enhancer target genes in top regions script
11.	EnrichR - genes with top UV ratio scores script
12.	EnrichR - enh target genes with top UV ratio scores script

Epigenetic feature analyses
1.	Epigenetic feature binning (1Mb bin)
2.	Histone mark & UV lesion correlation - FC and RankNorm script
3.	Euchr and Het scatterplot script
4.	PCA script

Chromatin states analyses
1.	Rank normalize UV lesion signal in chromatin state script
2.	UV lesion signal in chrom state boxplot script
3.	Histone mark binning to chromatin state script
4.	Histone mark heatmap script
5.	Genic feature and enhancer enrich in chrom state script
6.	Genic feature and enhancer heatmap script
7.	Repeat class enrich in chrom state script
8.	Repeat class heatmap
9.	Di-nucleotide frequency of each chromatin state script

## 3D Genome Modeling Script organization
Whole chromosome radial position scatterplot script
1. Whole chromosome radial pos scatterplot script

3D genome modeling scripts
1.	Chrom3D workflow batch script
2.	Test for optimal TAD overlap with top region script
3.	IMR TAD ID list for top regions script
4.	Batch script for coloring genome on Sherlock
5.	TAD nuclear distance boxplot script

## Melanoma Mutation Script organization

MELA-AU pre-processing scripts
1.	Subset MELA-AU dataset - only essential SSM info script
2.	Extracting donor IDs for cutaneous mela with UV sig script
3.	SSM dataset – cutaneous mela with UV sig only

Whole genome mutation rate scripts
1.	Genome melanoma mutation rate - 100kb and1Mb
2.	Mutation rate vs UV lesion abundance scatterplot script
3.	Mutation rate analysis of top susc. regions

Gene and enhancer mutation rate scripts
1.	Measurement of gene mutation rates
2.	Measurement of enhancer mutation rates
3.	Gene mutation rate analysis of top regions script
4.	Enhancer mutation rate analysis of top regions script

Cancer driver gene list scripts
1.	Dataset preparation script
2.	Heatmap – Cancer driver genes with top 64PP ratio scores
3.	Heatmap – Cancer driver genes with top CPD ratio scores
4.	Heatmap – Enhancer driver genes with top 64PP ratio scores
5.	Heatmap – Enhancer driver genes with top CPD ratio scores

## Nucleotide Excision Repair Script organization

Whole genome repair analyses scripts - 1Mb
1.	Binning XR-seq bw to 1Mb bin size
2.	Merging replicates
3.	Averaging plus and minus strand
4.	Master table with all time points 
5.	Cumulative repair vs UV lesion abundance scatterplot script

Chromatin state repair analysis scripts
1.	Binning XR-seq bw to chromatin state bin size
2.	Averaging plus and minus strand
3.	Master table with all time points
4.	Cumulative repair in chromatin state boxplots

Gene repair analysis scripts
1.	Subset uniquely mappable XRseq reads
2.	Measure XRseq signal in protein coding genes
3.	Averaging plus and minus strand
4.	Cumulative repair in protein coding genes 
5.	Gene repair analysis of top regions

Enhancer repair analysis scripts
1.	Measure unique XRseq signal in enhancer regions
2.	Average plus & minus strand repair signals
3.	Enhancer repair analysis of top susc. regions

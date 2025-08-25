setwd(r'{C:\Users\USER\Documents\BC_miRNA}')
dir.create("output")
dir.create("plots")

library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)
library(DESeq2)
library("pheatmap")
library("dichromat")
library("RColorBrewer")
library(ComplexHeatmap)
library(circlize)
library(stats)
library(openxlsx)
library(igraph)
library(Rgraphviz)
library(gridExtra)


#load data
data= read.csv("input/batch_corrected_data_2sets_18_unnormalized.csv") %>% column_to_rownames("X")

metadata= read.csv( "input/metadata_all.csv" )
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
metadata$Study_ID %>% unique()
metadata$Condition %>% unique()

condition= metadata$Condition[match(colnames(data), metadata$Sample_ID)]
condition= as.factor(condition)
table(condition)


#' ## #' ## construct DESeq dataset for disease status controlling for the source of sample/ biological batch effect:
dataset <- DESeqDataSetFromMatrix(countData = data,
                                  colData = metadata,
                                  design = ~ condition +Study_ID )

#' ## #' ## filter out genes with low reads
dataset <- dataset[ rowSums(counts(dataset)) > 500, ]
#' ## apply deseq model 
dds <- DESeq(dataset)

## get normalized count
### Deseq2 Normalization basis
### The size factor is calculated by taking the median of the ratios of observed counts to geometric means of counts for each gene. This means that for each sample,
### it measures how much the counts differ from the typical gene across all samples, which accounts for differences in sequencing depth.
### Once the size factors are calculated, the raw counts for each gene in each sample are divided by the size factor of that sample. 
#' ## get normalized count
norm.cts <- counts(dds, normalized=TRUE)

#write.csv(norm.cts, "output/all_rna_normalized_deseq.csv")
#' ## get coefficints of the model
coef(dds) %>% head()
#' ##  to bring the logFC of a gene to the general mean by estimating the overall dispersion of all genes in certain interval
res <- lfcShrink(dds, coef=  "condition_Tumor_vs_Normal" , type="apeglm")
#' ## summary for distribution of DE-genes
summary(res)

#' ## how many genes of adj p-values <= .05
sum(res$padj <= 0.05 & abs(res$log2FoldChange) > log(2, base=2), na.rm=TRUE)
sig= as.data.frame(res) %>% filter(., res$padj <= 0.05 & abs(res$log2FoldChange) >= log(2, base=2))

write.csv(as.data.frame(sig), "output/mrna_Deseq_res_fc2.csv")
write.csv(as.data.frame(res), "output/mrna_Deseq_res_all.csv")









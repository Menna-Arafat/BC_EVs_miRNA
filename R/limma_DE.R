setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

 load libraries
library(limma)
library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)
library("pheatmap")
library("dichromat")
library("RColorBrewer")
library(ComplexHeatmap)
library(circlize)

#load data
list.files(paste0(getwd(), "/output"))
# Load sample metadata and define condition/subtype groups.

study= "1"
# Set global variables
# 1st
if(study== "1"){
input= "/input1/"
suffix= "BC_vs_Ctrl"
plot_title= "Breast Cancer and Control Groups"
# expression data
data= read.csv( "output/data_all_combined_unnormalized_BC_vs_Ctrl.csv") %>% column_to_rownames("X")
# top features of sPls-DA MINT
topf= read.csv("output/BC_vs_Ctrl_top_features_component1.csv")

# metadata
metadata= read.csv( "metadata/metadata_all.csv" )
# Adjust metadata
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
unique(metadata$Study_ID)
#metadata levels
metadata$levels= metadata$Condition
lvl= c( "Ctrl", "BC")
table(metadata$levels)
metadata$levels <- factor(metadata$levels, levels= lvl)
metadata$Study_ID <- factor(metadata$Study_ID)
}else{
# 2nd
input= "/input2/"
sufffix= "subtypes"
plot_title= "Breast Cancer Subtypes"
data= read.csv("output/data_all_combined_unnormalized_subtypes.csv") %>% column_to_rownames("X")
# top features of sPls-DA MINT
topf= read.csv("output/subtypes_top_features_component1.csv")
# metadata
metadata= read.csv( "metadata/metadata_all.csv" )
# Adjust metadata
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
unique(metadata$Study_ID)
#metadata levels
metadata$levels= metadata$Subtype
lvl= unique(metadata$levels)
table(metadata$levels)
metadata$levels <- factor(metadata$levels, levels= lvl)
metadata$Study_ID <- factor(metadata$Study_ID)
}
#-------------------------------------------------------------------------------
# Build the design matrix
design <- model.matrix(~ 0+ Study_ID +levels,
                       data = metadata,
                       contrasts.arg = list(levels = contrasts(metadata$levels , contrasts = FALSE))) %>% as.data.frame()
head(design)
colnames(design)= gsub("levels|Study_ID","", colnames(design))
#-------------------------------------------------------------------------------
#Run Limma Pipeline
#One subtype vs all others, controlling for other covariates as study id
for(i in  unique(metadata$levels) ){
  des= design[, c( as.character(unique(metadata$Study_ID)), i)]
  #removes columns where all values are 0 (i.e., no samples in the group).
  des <- des[, colSums(des) > 0]

  v <- voom(data, des)
  
  # Fit the linear model
  fit <- lmFit(v, des)
  # Apply empirical Bayes moderation with trend adjustment
  fit <- eBayes(fit, trend = TRUE)
 
    Get the top differentially expressed genes
  topTable <- topTable(fit, coef = i,number = Inf)
  sig= topTable[topTable$P.Value < .05 & abs(topTable$logFC) > log2(1.5) ,]
  
  write.csv(sig, paste0("output/DE_limma_sig_" , i, ".csv"))
  write.csv(topTable , paste0("output/DE_limma_all_" , i, ".csv"))
  
  #Intersect top features sPLSDA with DE
  y= intersect(row.names(sig), topf$X)
  write.csv(y,  paste0("output/intersected_sig_" , i, ".csv"), row.names = F)
  
}
#-------------------------------------------------------------------------------
# #Get the results for a specific contrast
# coef(fit) %>% head()
# contrast <- makeContrasts(Tumor - Normal, levels = design)
# fit2 <- contrasts.fit(fit, contrast)
# fit2 <- eBayes(fit2)
# 
#  # Get the top differentially expressed genes
# topTable <- topTable(fit2, number = Inf)

  
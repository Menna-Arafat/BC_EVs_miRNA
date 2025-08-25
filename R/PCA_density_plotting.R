setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

library(graphlayouts)
library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)
library(igraph)
library(Matrix)
library(ggplot2)
library(plotly)
library(plotly)
library(readxl)
library(ggpubr)
#remotes::install_github("EvaYiwenWang/PLSDAbatch")
library(PLSDAbatch)
library(mixOmics)

# get data
list.files(paste0(getwd(), "/output"), pattern= "batch_corrected")
files= c("batch_corrected_data_BC_vs_Ctrl.csv" ,"batch_corrected_data_subtypes.csv"  )
i=1
data= read.csv(paste0(getwd(), "/output/", files[i])) %>% column_to_rownames("X")
metadata= read.csv( "metadata/metadata_all.csv" )
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
metadata$Study_ID %>% unique()
metadata$levels= metadata$Condition
group_colors = c( "#FFA500CC", "darkslateblue")
names(group_colors)= unique(metadata$levels)
plot_title= "Breast Cancer and Control Groups"
suffix="BC_vs_Ctrl"


# Perform PCA
# scale data 
data= as.data.frame(t(scale(t(data), center = TRUE, scale = TRUE )))
pca_res <- pca(t(data), ncomp = 3, center = FALSE, scale = FALSE)

# Get sample coordinates
pca_coords = pca_res$x[, 1:2]  # Take only PC1 and PC2
rownames(pca_coords) = colnames(data)

# Get variance explained
expl.var = summary(pca_res)$importance[2, 1:2] * 100  # PC1 and PC2 %

# Define batch and treatment variables
batch = metadata$Study_ID[match(rownames(pca_coords), metadata$Sample_ID)]
trt = metadata$levels[match(rownames(pca_coords), metadata$Sample_ID)]

names(batch)= names(trt)= colnames(data)
# Call the custom plotting function
Scatter_Density(
  pca_res,
  batch = batch,
  trt = trt,
  xlim = c(-4.5, 5),
  ylim = c(-3, 4),
  batch.legend.title = 'Study',
  trt.legend.title = 'Group',
  title = paste0("PCA with Density: ", plot_title)
)

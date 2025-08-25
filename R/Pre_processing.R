
setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

#load libraries
library(limma)
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
library(magick)

list.files(paste0(getwd(), "/output"))
# set global variables

study= "1"

if(study== "1"){
# 1st
input= "/input1/"
suffix= "BC_vs_Ctrl"
plot_title= "Breast Cancer and Control Groups"
}else{
# 2nd
input= "/input2/"
suffix= "subtypes"
plot_title= "Breast Cancer Subtypes"
}

#-------------------------------------------------------------------------------
#load expression data
#list.files(paste0(getwd(), "/input"))
list=list()
for (i in list.files(paste0(getwd(), input), pattern= ".csv")){
  list[[i]]= read.csv(paste0(getwd(), input, i), sep= ";") #, sep= ";"
  print(dim(list[[i]]))
  
}

#list[[1]]= read.csv(paste0(getwd(), input, "GSE141326_Raw_count_with_ID_names.csv"))
lapply(list, function(x) colnames(x))

# remove extra colmns
list= lapply(list, function(x) x= x[, !colnames(x) %in% 
                                      c( "Accession..", "Class.Name", "Analyte.Type")])

#excluded for low input miRNA
#x2= read.delim("input/GSE222681_Raw_Count.txt" )
#x3= read.delim( "input/GSE239341_Raw_Count.txt" )
# x4= read.delim(  "input/GSE270497_All.Counts.exp.txt") %>% dplyr::select(-c(1,3)) 
# #adjust mirna names
# x4$smallRNApreName = gsub("mir", "miR", x4$smallRNApreName)
# x4= x4 %>% group_by(smallRNApreName)%>% 
#   slice_max(rowSums(across(where(is.numeric))), n = 1, with_ties = FALSE) %>%
#   ungroup()
#x5= read.delim("input/GSE256523_Normalised_EV_count_file.tabular.txt")
#-------------------------------------------------------------------------------
# preprocessing: Standardize miRNA IDs && Remove duplicates (keep max expressed).

list= lapply(list, function(x){ 
  colnames(x)[1]= "miRNA"
  # retain first 3 parts of ID
  x$miRNA= gsub("^((\\w+-){2}\\w+).*", "\\1", x$miRNA)
  # change data type to numeric
  x[,-1]= as.data.frame(lapply(x[,-1], as.numeric))
  # remove duplicates by keeping the highest expression of both
  x = x %>% group_by(miRNA)%>% 
    slice_max(rowSums(across(where(is.numeric)), na.rm=TRUE), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    column_to_rownames("miRNA")
  
  print(dim(x))
  return(x)
})
#-------------------------------------------------------------------------------
# Normalization
# list_z <- lapply(list, function(x) {
#   shared <- Reduce(intersect,  lapply(list, rownames))
#   x= x[row.names(x) %in% shared, ]
#   x= as.data.frame(t(scale(t(x), center = TRUE, scale = TRUE )))
#   return(x)
# })
# # Combine all 
# data <- do.call(cbind, list_z)

# Without normalization
# Data Combination
# Keep shared genes only in all datasets
shared_genes= Reduce(intersect, lapply(list, rownames))

#Subset each dataset by shared genes
list= lapply(list, function(df) df[shared_genes, , drop = FALSE])

# Combine all 
names(list)= NULL
data <- do.call(cbind, list)
dim(data)

# Filteration 
# Filter genes with too many NAs (>60% missing)
threshold <- 0.6
na_fraction <- rowMeans(is.na(data))
data <- data[na_fraction <= threshold, ]
# Impute remaining Nas with 0
data[is.na(data)]=0
dim(data)

# export
write.csv(data, paste0("output/data_all_combined_unnormalized_", suffix, ".csv"))
#-------------------------------------------------------------------------------
#data= read.csv("output/data_all_combined_unnormalized_BC_vs_Ctrl.csv") %>% column_to_rownames("X")
#data= read.csv("input/data_all_combined_unnormalized_subtypes.csv") %>% column_to_rownames("X")

# Batch effect Correction
# 1st
if(study== "1"){
    metadata= read.csv( "metadata/metadata_all.csv" )
    metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
    setdiff(colnames(data), metadata$Sample_ID )
    metadata= metadata[metadata$Sample_ID %in% colnames(data),]
    metadata$Study_ID %>% unique()
    metadata$levels= metadata$Condition
    group_colors = c( "#FFA500CC", "darkslateblue")
    names(group_colors)= unique(metadata$levels)
}else{
#2nd
    metadata= read.csv( "metadata/metadata_all.csv" )
    metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
    setdiff(colnames(data), metadata$Sample_ID )
    metadata= metadata[metadata$Sample_ID %in% colnames(data),]
    metadata$Study_ID %>% unique()
    #metadata levels
    metadata$levels= metadata$Subtype
    table(metadata$levels)
    set.seed(1234)
    group_colors <- colorRampPalette(c( "#180C3CFF", "#961915", "maroon3" ,"purple",  "#01665E" ,
                                        "#557192" , "#B0C4DE", "#E3B31C", "grey32" ))(length(unique(metadata$levels))) 
    
    names(group_colors)= unique(metadata$levels)
}

# Apply limma to remove the batch effect of the study factor 
# and export residuals as our biological signal
batch <- metadata$Study_ID[match(colnames(data), metadata$Sample_ID)]
#unique(batch)

# Design matrix for biological conditions 
design <- model.matrix(~1, data = metadata)
head(design) # just intercept
# Apply removeBatchEffect from limma
# removeBatchEffect subtracts the batch effect from the original data, 
#resulting in a residual matrix that contains condition effect. 
residuals_matrix <- removeBatchEffect(data, batch = batch, design = design)

# Export residuals (biological signal)
write.csv(residuals_matrix, file = paste0("output/batch_corrected_data_", suffix, ".csv"))
#-------------------------------------------------------------------------------
# PCA
# Function
plot.pca= function(data, metadata, plot_title, suffix){
# scale data 
data= as.data.frame(t(scale(t(data), center = TRUE, scale = TRUE )))

# Reduce dimensionality using PCA
pca_res= prcomp(t(data), center =FALSE, scale. = FALSE) #samples in rows

#pca scores for samples
PCA_data = pca_res$x %>% as.data.frame()

#variance explained by PCs
vars <- apply(PCA_data, 2, var)
PC1_var <- round( ( var(PCA_data$PC1) / sum(vars) ) * 100, 1) 
PC2_var <- round( ( var(PCA_data$PC2) / sum(vars) ) * 100, 1)

# Variance explained by each principal component
explained_variance <- summary(pca_res)$importance[2, ]
head(explained_variance)

# Add group information
PCA_data$Condition <-  metadata$levels[match(colnames(data), metadata$Sample_ID)]
head(PCA_data$Condition)
table(PCA_data$Condition)
# Add source of each dataset
PCA_data$Study= metadata$Study_ID[match(colnames(data), metadata$Sample_ID)]
#table(PCA_data$Study)

p= ggplot(PCA_data, aes(x = PC1, y = PC2,  shape = Study , color = Condition)) + 
  geom_point(size = 2) +  
  scale_color_manual(values =group_colors, name = "Groups") +
  ggtitle(paste0("PCA Plot of ", plot_title))  +
  #ggtitle("Breast Cancer Subtype Stratification in PCA Space") +
  xlab(paste0("PC1 ", "(", PC1_var , "%" ,")")) +
  ylab(paste0("PC2 ", "(", PC2_var , "%", ")")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    #panel.border = element_blank()
  ) 

print(p)
ggsave(paste0("plots/PCA_plot_", suffix, ".png"), p, width = 9, height = 7)
}

#apply
data= residuals_matrix 
plot.pca(data, metadata, plot_title, paste0(suffix, ".after"))
#===============================================================================
# Convert ensembl ids to miRNA mirbase ID
#convert ids
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listFilters(ensembl)      # To find filter names like 'ensembl_gene_id'
listAttributes(ensembl)   # To find attributes like 'mirbase_id'
# Example vector of Ensembl gene IDs
ensembl_ids = x5$GeneID %>% unique()

# Get mapping
miRNA_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "mirbase_id"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# match id to ensmble
x5$ID= miRNA_map$external_gene_name[match(x5$GeneID, miRNA_map$ensembl_gene_id) ]
#-------------------------------------------------------------------------------
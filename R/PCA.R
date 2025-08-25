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
    geom_point(size = 2)+
    scale_color_manual(values = group_colors, name = "Groups") +  
    ggtitle(paste0("PCA Plot of ", plot_title))  +
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

# get data
list.files(paste0(getwd(), "/output"), pattern= "batch_corrected")
files= c("batch_corrected_data_BC_vs_Ctrl.csv" ,"batch_corrected_data_subtypes.csv"  )

#apply
for(i in 1:length(files)){
  if(i == 1){
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
  }else{
    #2nd
    data= read.csv(paste0(getwd(), "/output/", files[i])) %>% column_to_rownames("X")
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
    plot_title= "Breast Cancer Subtypes"
    suffix= "Subtypes_"
  }
  
plot.pca(data, metadata, plot_title, paste0(suffix, ".after"))



}

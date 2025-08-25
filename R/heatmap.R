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
data= read.csv("input/batch_corrected_data_3sets_subtyping.csv") %>% column_to_rownames("X")

DE_list= list()
for(i in list.files(paste0(getwd(), "/output/sig"))){
  DE_list[i] = read.csv(paste0("output/sig/", i))
  }
 
DE= DE_list %>% unlist() %>% unname() %>% unique() 
topf= read.csv("output/top_features_component1_subtypes.csv")

y= intersect(DE, topf$X)
write.csv(y, "output/intersected_all_DE_subtypes_splsda.csv", row.names = F)
# DE= read.csv("output/DE_limma_sig_BC_vs_Ctrl.csv")
# DE= DE %>% arrange(desc(abs(logFC))) %>% slice(1:40)
heat_data= data[ match(y, rownames(data)) ,]
heat_data= scale(heat_data)

metadata= read.csv( "input/metadata_all.csv" )
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
metadata$Study_ID %>% unique()
metadata$Condition %>% unique()

metadata$Condition <- factor(metadata$Condition, levels = c("BC", "Ctrl"))
metadata$Subtype <- factor(metadata$Subtype, levels = unique(metadata$Subtype))

table(metadata$Condition)
table(metadata$Subtype)


#plotting
unique( metadata[, c("Condition")])

ta <- HeatmapAnnotation(
  Group = metadata$Subtype,
  col = list(
    #Group = c("BC" = "#CB2892", "Ctrl" = "#B7DD34")
    Group = c("Luminal A" =  "#BF62AC","Luminal B"= "blue", "TNBC"= "red", 
              "Her2"=  "#007E81FF" , "Luminal HER2"= "#FDE725", "Ctrl" = "#B7DD34")
  ),
  annotation_height = unit(10, "mm"),
  annotation_legend_param = list(
    title = "Group",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 10)
  )
)




palt1= colorRampPalette(c( "lightyellow2" ,"#DFC27D" , "#FCAA0FFF", "darkred"))(256)
palt2= colorRamp2(c(-2, 0, 2), c( "#483D8B","lightyellow2",  "orange2"))


heatmap <- Heatmap(
  matrix = as.matrix(heat_data),
  name = "Z scores",
  col = palt2,
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  top_annotation = ta,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),   # Adjusts the title font size
    labels_gp = gpar(fontsize = 10),  # Adjusts the label font size
    grid_width = unit(9, "mm"),       # Increases the width of legend color blocks
    grid_height = unit(9, "mm")       # Increases the height of legend color blocks
  )
)

heatmap 
#' ## Save 
png("plots/heatmap_subtypes_DE_intersected.png",width = 6000, height = 7000, res = 600)
draw(heatmap, annotation_legend_side =  "right")
dev.off()


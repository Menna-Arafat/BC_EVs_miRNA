# Tutorial https://mixomics.org/mixmint/mint-stem-cells-case-study-2/
setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

library(graphlayouts)
library(igraph)
library(WGCNA)
library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)
library(igraph)
library(Matrix)
library(ggplot2)
library(plotly)
library(plotly)
library("mixOmics") 
library(readxl)
library(ggpubr)

list.files(paste0(getwd(), "/output"))
# Set global variables
prefix= "BC_vs_Ctrl_"
plot_title= "Breast Cancer and Control Groups"

# Adjust metadata
data= read.csv("output/batch_corrected_data_BC_vs_Ctrl.csv") %>% column_to_rownames("X")

metadata= read.csv( "metadata/metadata_all.csv" )
metadata$Sample_ID= gsub("Positive\\.|Negative\\.", "",metadata$Sample_ID )
setdiff(colnames(data), metadata$Sample_ID )
metadata= metadata[metadata$Sample_ID %in% colnames(data),]
metadata$Study_ID %>% unique()
metadata$Condition %>% unique()

condition= metadata$Condition[match(colnames(data), metadata$Sample_ID)]
condition= as.factor(condition)
table(condition)

#BC Ctrl 
#60   42 
study= metadata$Study_ID[match(colnames(data), metadata$Sample_ID)]
study= as.factor(study)
table(study)

# PLS-DA pipeline
data.t= t(data) # samples in row
dim(data.t) 

# nearZeroVar function (should be set to TRUE in particular for data with many zero values)
#keepX to be a vector of length equal to ncomp, i.e., one value per component.
ncomp =4

splsda.model <- mint.splsda(data.t, condition, study = study, ncomp =ncomp ,
                            near.zero.var=TRUE, keepX=rep(150,ncomp  )) 

#which features were near-zero variance, and got removed
splsda.model$nzv
ex= splsda.model$nzv$Metrics %>% row.names()
#data.t= data.t[, !colnames(data.t)%in% ex ]

#-------------------------------------------------------------------------------
# TUNING
# Number of componenets to keep
splsda.perf <- perf(splsda.model, near.zero.var=TRUE,  cpus = 3, validation = "Mfold", folds = 5, nrepeat = 10) # undergo performance optimisation

png(paste0("plots/", prefix, "classification_error.png"), width = 2000, height = 1600, res = 300)
plot(splsda.perf)
dev.off()

# Number of features to keep
# sPLS-DA tuned with optimal parameters
splsda.tune <- tune(data.t, condition, study = study, 
                    ncomp =ncomp , 
                    near.zero.var=TRUE,
                    test.keepX =rep(150,ncomp), 
                    method = 'mint.splsda', 
                    measure = 'overall', # balanced error rate
                    dist = "max.dist")

#Diamonds represent the optimal number of features on a given component. 
#Balanced error rate found on the vertical axis and is the metric to be minimised.
plot(splsda.tune, sd = FALSE)

optimal.keepX <- splsda.tune$choice.keepX # extract optimal values
optimal.keepX  # 150
#--------------------------------------------------------------------------------
# Final model with optimal number of features
#generate optimal model using tuned parameters
final.splsda.model <- mint.splsda(data.t, condition, study = study, 
                                  ncomp = 2, 
                                  keepX = 150)

# sPLS-DA plot
png(paste0("plots/", prefix, "splsda_global_plot.png"), width = 5000, height = 3600, res = 600)
plotIndiv(final.splsda.model, study = 'global', 
          legend = TRUE,
          title = paste0("MINT sPLS-DA for " , plot_title),
          subtitle = 'Global', ellipse = T)

dev.off()


png(paste0("plots/", prefix, "splsda_for_each_study.png"), width = 5000, height = 4600, res = 600)
plotIndiv(final.splsda.model, study = 'all.partial', 
          legend = TRUE,
          title =  paste0("MINT sPLS-DA for " , plot_title))
dev.off()
#-------------------------------------------------------------------------------
# heatmap for features
row_colors= color.mixo(as.numeric( condition))
unique(row_colors)
unique(condition)

png(paste0("plots/", prefix, "heatmap.png"), width = 4600, height = 4800, res = 600)
cim(final.splsda.model, comp = 1, margins = c(6, 9),
    row.names = FALSE,
    row.sideColors = row_colors,
    title = "MINT sPLS-DA, Comp. 1",
    legend = list(x = "right", legend = unique(condition), fill = unique(row_colors) ))
dev.off()

#png("plots/network_plot.png", width = 4600, height = 4800, res = 600)
# network(final.splsda.model, comp = 1, 
#         color.node = c(color.mixo(1), color.mixo(2)), 
#         shape.node = c("rectangle", "circle"))
# dev.off()
#-------------------------------------------------------------------------------
#generate optimal model using tuned parameters; optimal.keepX
final.splsda.model2 <- mint.splsda(data.t, condition, study = study, 
                                   ncomp = 2, 
                                   keepX = 150)

#Performance of the model
png(paste0("plots/", prefix, "ROC_plot.png"), width = 4600, height = 4800, res = 600)
auroc(final.splsda.model2, roc.comp = 1, print = FALSE, plot = TRUE)
dev.off()



top_features <- selectVar(final.splsda.model, comp = 1)
top_features= top_features[["value"]]
#export
write.csv(top_features, paste0("output/", prefix, "top_features_component1.csv"))


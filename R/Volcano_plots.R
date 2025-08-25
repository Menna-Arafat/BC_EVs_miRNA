setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("plots/volcanos/")
#' ## BiocManager::install("fields")
library(tibble)
library("RColorBrewer")
library("circlize")
library(biomaRt)
library(dplyr)
library(plyr)
library(fields)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggrepel)

plot.volcano <- function(output_dir = "output", plot_dir = "plots/volcanos", fc_threshold = 1.5) {

  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # List DE result files
  files <- list.files(output_dir, pattern = "DE_limma_all", full.names = TRUE)
  
  for (file_path in files) {
    # Extract the name (group/subtype) from filename
    name <- gsub(".*_([[:alnum:] ]+)\\.csv$", "\\1", basename(file_path))
    
    # Read the DE result
    res <- read.csv(file_path)
    res <- as.data.frame(res)
    
    # Define direction
    res$Direction <- ifelse(res$P.Value <= 0.05 & res$logFC >= log2(fc_threshold), "Up",
                            ifelse(res$P.Value <= 0.05 & res$logFC <= -log2(fc_threshold), "Down",
                                   "Non-significant"))
    
    # Keep only labels for up or down expressed features
    res <- res %>% mutate(features = ifelse(Direction %in% c("Up", "Down"), X, ""))
    
    # Compute -log10(p-value)
    res$log10p <- -log10(res$P.Value)
    
    # Plot parameters
    xminma <- -3.5
    xmaxma <- 3.5
    
    # Generate the plot
    volcano <- ggplot(res, aes(x = logFC, y = log10p, color = Direction, label = features)) +
      geom_point(size = 1.2, alpha = 0.7) +
      geom_rug(alpha = 0.6) +
      scale_color_manual(values = c("Up" = "red2", "Down" = "darkslateblue", "Non-significant" = "grey66")) +
      xlab('log2 Fold Change') +
      ylab('-log10 p-value') +
      scale_x_continuous(limits = c(xminma, xmaxma)) +
      theme_bw() +
      theme(legend.title = element_blank()) +
      geom_vline(xintercept = c(-log2(fc_threshold), 0, log2(fc_threshold)), linetype = c("dotted", "solid", "dotted")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_text_repel(aes(label = features),
                      size = 2,
                      max.overlaps = 10,
                      segment.color = "grey50",
                      color = "black") +
      ggtitle(paste0("Volcano Plot of DE miRNAs: ", name))
    
    # Print and save
    print(volcano)
    ggsave(filename = paste0(plot_dir, "/Volcanoplot_", name, ".jpg"),
           plot = volcano,
           dpi = 600,
           width = 7,
           height = 4)
  }
}

plot.volcano()

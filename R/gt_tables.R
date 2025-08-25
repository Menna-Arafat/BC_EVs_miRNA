
setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

#load libraries
library(gridExtra)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(scales)
library(gplots)
library(ggrepel)
library(ggpubr)
library(tibble)
library(cowplot)
library(gtsummary)
library("RColorBrewer")
library(kableExtra)
library(gt)
library(gtsummary)

# load data
table= read.csv("output/enrichment/enrichment__net_targets.csv")
table$ID= gsub("%.*", "", table$ID)
table$ID= substr(table$ID, 1, 60)

enrich_list=list()
for (i in unique(table$.id )){

  enrich_table= table[table$.id== i, ]
  enrich_x= enrich_table[enrich_table$p.adjust <= .05,] # 
  
  
  if(nrow(enrich_x)< 5){
    enrich_x= enrich_table[enrich_table$p.adjust <= .5,]
    
    if (nrow(enrich_x) > 100){
      enrich_x= enrich_x[1:100,]
    }
  }
  
  enrich_list[[i]]= enrich_x$ID
}
#-------------------------------------------------------------------------------
# Extract core pathways
tnbc_core_indices <- c(
  1,   # "C-MYB TRANSCRIPTION FACTOR NETWORK"
  2,   # "NEGATIVE REGULATION OF AUTOPHAGY"
  7,   # "G1/S TRANSITION OF MITOTIC CELL CYCLE"
  10,  # "PID_CMYB_PATHWAY"
  14,  # "INTEGRATED CANCER PATHWAY"
  36,  # "DNA DAMAGE RESPONSE"
  50,  # "TP53 NETWORK"
  55,  # "PI3K AKT SIGNALING"
  83   # "HIPPO SIGNALING"
)
her2_core_indices <- c(
  1,   # "G1 TO S CELL CYCLE CONTROL"
  3,   # "PID_RB_1PATHWAY"
  5,   # "MIRNA REGULATION OF DNA DAMAGE RESPONSE"
  9,   # "RETINOBLASTOMA GENE IN CANCER"
  11,  # "ONCOGENE INDUCED SENESCENCE"
  24,  # "PID_MYC_ACTIV_PATHWAY"
  30,  # "OXIDATIVE STRESS INDUCED SENESCENCE"
  41,  # "PID_P53_DOWNSTREAM_PATHWAY"
  54,  # "P53 PATHWAY"
  68   # "PID_BETA_CATENIN_NUC_PATHWAY"
)
luminalA_core_indices <- c(
  2,   # "CELLULAR RESPONSE TO HEAT STRESS"
  6,   # "REGULATION OF HSF1-MEDIATED HEAT SHOCK RESPONSE"
  30,  # "PI3K AKT SIGNALING"
  99   # "CHROMATIN ORGANIZATION"
)
luminalHer2_core_indices <- c(
  1,   # "PROTEIN PHOSPHORYLATION"
  7,   # "ESR-MEDIATED SIGNALING"
  14,  # "TRANSCRIPTIONAL REGULATION BY TP53"
  17,  # "TP53 REGULATES METABOLIC GENES"
  20,  # "CLASS I PI3K SIGNALING EVENTS MEDIATED BY AKT"
  24,  # "RETINOBLASTOMA GENE IN CANCER"
  29,  # "G1 TO S CELL CYCLE CONTROL"
  53,  # "ESTROGEN-DEPENDENT GENE EXPRESSION"
  92   # "C-MYB TRANSCRIPTION FACTOR NETWORK"
)

core_her2_pathways <- enrich_list[["Her2"]][her2_core_indices]
core_luminalA_pathways <- enrich_list[["Luminal A"]][luminalA_core_indices]
core_luminalB_pathways <- enrich_list[["Luminal B"]]
core_luminalher2_pathways <- enrich_list[["Luminal Her2"]][luminalHer2_core_indices]
core_TNBC_pathways <- enrich_list[["TNBC"]][tnbc_core_indices]


results_table1= table[table$ID %in% c(core_her2_pathways,core_luminalA_pathways,
                                     core_luminalher2_pathways,core_TNBC_pathways) &
                                     table$p.adjust <=.05, ]

results_table2= table[table$ID %in% core_luminalB_pathways & table$.id== "Luminal B",]
results_table= rbind(results_table1, results_table2)
results_table= results_table %>% slice(-c(1,2)) %>% dplyr::select(1:2)
head(results_table)

res.tbl=  results_table %>%
      group_by(.id) %>%
      mutate(row = row_number()) %>%
      pivot_wider(names_from = .id, values_from = ID) %>%
      select(-row)

res.tbl[is.na(res.tbl)]=""
head(res.tbl)

#visualize
library(gt)

tbl <- gt(res.tbl) %>%
  tab_header(
    title = md("**Core Pathways across Breast Cancer molecular subtypes**"),
    subtitle = md("Over-representation Enrichment Analysis")
  ) %>%
  # Style the title
  tab_style(
    style = cell_text(weight = "bold", size = "large"),
    locations = cells_title(groups = "title")
  ) %>%
  # Style the subtitle
  tab_style(
    style = cell_text(size = "medium", color = "black"),
    locations = cells_title(groups = "subtitle")
  ) %>%
  # Style the body
  tab_style(
    style = cell_text(style = "italic", color = "darkblue"),
    locations = cells_body()
  ) %>%
  opt_table_font(
    font = list(
      google_font("IBM Plex Sans"),
      default_fonts()
    )
  ) %>%
  opt_row_striping() %>%
  # Change Column name
  cols_label(
     Her2 = "HER2-Enriched",
    `Luminal A` = "Luminal A",
    `Luminal HER2` = "Luminal HER2",
    TNBC = "Triple-Negative (TNBC)",
    `Luminal B` = "Luminal B"
  ) %>%
  tab_options(
    table.border.top.color = "gray80",
    table.border.bottom.color = "gray80",
    row.striping.background_color = "#f9f9f9",
    heading.align = "left",
    column_labels.font.weight = "bold",
    table.font.size = 14
  ) %>%
  tab_footnote(
    footnote = "Adj.P.Value ≤ 0.05 (Subtypes), P.Value ≤ 0.05 (Luminal B)"
  ) %>%
  #  Style the footnote
  tab_style(
    style = cell_text(size = "small", weight = "bold"),
    locations = cells_footnotes()
  )


# Print the table
print(tbl)
gtsave(tbl, "plots/tbl_pathways.png")

#///////////////////////////////////////////////////////////////////////////////
# BC vs Ctrl
sig.p=enrich_list[["BC"]]
res.tbl= table[table$ID %in% sig.p & table$.id== "BC",]
res.tbl$ID
res.tbl=res.tbl %>% dplyr::select(ID, geneID) %>% slice(-c(1,2,6, 7,9,14, 15,16,22,29, 32, 33))


tbl2 <- gt(res.tbl) %>%
  tab_header(
    title = md("**Core Biological Processes in Breast Cancer**"),
    subtitle = md("Over-representation Enrichment Analysis")
  ) %>%
  # Style the title
  tab_style(
    style = cell_text(weight = "bold", size = "large"),
    locations = cells_title(groups = "title")
  ) %>%
  # Style the subtitle
  tab_style(
    style = cell_text(size = "medium", color = "black"),
    locations = cells_title(groups = "subtitle")
  ) %>%
  # Style the body
  tab_style(
    style = cell_text(style = "italic", color = "darkblue"),
    locations = cells_body()
  ) %>%
  opt_table_font(
    font = list(
      google_font("IBM Plex Sans"),
      default_fonts()
    )
  ) %>%
  opt_row_striping() %>%
  # Change Column name
  cols_label(
    ID = "Biological Processes",
    geneID= "Gene.ID"
  ) %>%
  tab_options(
    table.border.top.color = "gray80",
    table.border.bottom.color = "gray80",
    row.striping.background_color = "#f9f9f9",
    heading.align = "left",
    column_labels.font.weight = "bold",
    table.font.size = 14
  ) %>%
  tab_footnote(
    footnote = "Adj.P.Value ≤ 0.05"
  ) %>%
  #  Style the footnote
  tab_style(
    style = cell_text(size = "small", weight = "bold"),
    locations = cells_footnotes()
  )


# Print the table
print(tbl2)
gtsave(tbl2, "plots/tbl_pathways_BC.png")

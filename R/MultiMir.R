# Mirna targets

setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

#load libraries
#BiocManager::install("multiMiR")
library(multiMiR)
library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)

# Find miRNA targets using Multimir from Searching mirecords, mirtarbase, tarbase databased
# intiate output list
results= list()

for (i in list.files(paste0(getwd(), "/output"), pattern = "intersected") ){
  
  name= gsub(".*_([[:alnum:] ]+)\\.csv$", "\\1", i)
  x= read.csv(paste0("output/", i)) %>% dplyr::select("x") %>% 
                                        unlist() %>% unname()

  #run multimir
  multimir_results= get_multimir( org     = 'hsa',
                                  mirna   = x,
                                  table   = 'validated',
                                  summary = TRUE)
  
  results[[name]] <- multimir_results@data[3:4]
  
}

res_df= do.call(rbind, results)
res_df= res_df %>%  rownames_to_column("Group")
res_df$Group=  gsub("\\.\\d+", "", res_df$Group)
head(res_df)
write.csv(res_df, "output/mir_targets_de_intersect_subtypes.csv", row.names = F)





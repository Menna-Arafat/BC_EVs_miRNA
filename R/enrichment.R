setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")

library(stats)
library(tibble)
library(plyr)
library(dplyr)
library(tidyverse)
library(fgsea)
library("GSEABase")
library("EnrichmentBrowser")
library(clusterProfiler)
library(stats)
library(tibble)
library("RColorBrewer")
library("circlize")
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

# gmt files
gmt_mirna= gmtPathways("c:/Users/USER/Documents/resources/TAM_Function_precursor_miEAA 2.0_resources.gmt")
gmt_all = gmtPathways("c:/Users/USER/Documents/resources/Human_GOBP_AllPathways_noPFOCR_with_GO_iea_July_01_2024_symbol.gmt.txt")

#load data
#[1]
de_list=list()
for (i in list.files(paste0(getwd(), "/output/DE/"), pattern = "limma_sig") ){
  de_list[i]= read.csv(paste0("output/DE/", i))
}

names(de_list)= gsub(".*_([[:alnum:] ]+)\\.csv$", "\\1",names(de_list))


#[2]
# get overlapped targets that are tageted by at least 3 differentially expressed miRNAs
targets= read.csv("output/mir_targets_de_intersect_subtypes.csv")
unique(targets$Group)
head(targets)
list <- split(targets[, -1], targets$Group)

target_list=list()
for(i in names(list)){
  net = list[[i]]
  colnames(net)= c("ID", "Target")
  net= net[net$ID != ""& net$Target != "", ]
  
  if(i %in% c( "Luminal B", "Luminal A"  )){
    num=1
  }else{
    num=3
  }
  # get overlapped targets that are tageted by at least 3 differentially expressed miRNAs
  shared_targets <- net %>%
                    group_by(Target) %>%
                    summarise(n_miRNAs = n_distinct(ID)) %>%
                    filter(n_miRNAs >= num)
  
  #subset
  net <- net %>%
         filter(Target %in% shared_targets$Target)
  
  genes = net %>% unlist() %>% unname() %>% unique()
  target_list[[i]]= genes
}
#-------------------------------------------------------------------------------
data= "1"

if(data== "1"){
gmt= c(gmt_mirna)
en_list= de_list
name_plot= "de_mirna"

}else{
gmt= c(gmt_mirna, gmt_all)
en_list= target_list
name_plot= "_net_targets"

}


#filter gene set to have at least 5 terms pr set
gset_sub= gmt[lapply(gmt, length) >= 2& !grepl("^PMC",names(gmt)) ]
# convert list to long formats
gmt_long= stack(gset_sub) 
gmt_long= gmt_long[,c(2,1)]  
gmt_long$values= gsub("^((\\w+-){2}\\w+).*", "\\1", gmt_long$values)
gmt_long$values= gsub("mir", "miR", gmt_long$values)
head(gmt_long)
tail(gmt_long)


# Run clusterprofiler
#Pathway enrichment clusterprofiler
enrich=  lapply(en_list, function(i){
  #GSEA
  #i = i[!duplicated(names(i))]  
  #e = GSEA(i, TERM2GENE = gmt_long, pvalueCutoff = 1, minGSSize = 5)
  #ORA
  e= enricher(i, TERM2GENE =gmt_long, minGSSize = 1, pvalueCutoff = .5,  qvalueCutoff = .5)
  return(as.data.frame(e))
})

sig= lapply(enrich, function(i){
  i= i[i$p.adjust <= .2, ]
  return(as.data.frame(i))
})

sig_df=ldply(sig, rbind)

write.csv(sig_df, paste0("output/enrichment_", name_plot ,".csv"), row.names = F)
#-------------------------------------------------------------------------------

data= "2"
# dot plot
if(data== "1"){
  dir.create("plots/enrich_DE")
  name2= "DE_"
  dir= "plots/enrich_DE/"
  table= read.csv("output/enrichment_de.csv")

}else{
  dir.create("plots/enrich_targets")
  name2= "net_targets_" 
  dir= "plots/enrich_targets/"
  table= read.csv("output/enrichment__net_targets.csv")
}

table$.id= as.factor(table$.id)
unique(table$.id)

# Dot plot
for (i in unique(table$.id )){
      name= gsub(".*_([[:alnum:] ]+)\\.csv$", "\\1", i)
      enrich_table= table[table$.id== i, ]
      enrich_x= enrich_table[enrich_table$p.adjust <= .05,]
        
      
      if(nrow(enrich_x)< 5){
        enrich_x= enrich_table[enrich_table$p.adjust <= .5,]
      }
      
      if (nrow(enrich_x) > 40){
        enrich_x= enrich_table[1:40,]
      }
      
      enrich_table= enrich_x
      #fraction read as factor
      enrich_table$GeneRatio= sapply(enrich_table$GeneRatio, function(x) eval(parse(text=x)))
      enrich_table$fold_enrichment= enrich_table$GeneRatio
      enrich_table$fold_enrichment
      #enrich_table$fold_enrichment = enrich_table$intersection_size/ enrich_table$term_size
      colnames(enrich_table)[2] <- "new_pathway_name"
      enrich_table$new_pathway_name= gsub("%.*", "", enrich_table$new_pathway_name)
      enrich_table$new_pathway_name= substr(enrich_table$new_pathway_name, 1, 50)
      enrich_table$P. <- -log10(enrich_table$p.adjust) %>% round(2)
      enrich_table$P.
      
      plot1 = ggplot(data=enrich_table, aes( y = reorder(new_pathway_name, P.), 
                                             x = P., 
                                             size = fold_enrichment, 
                                             color = p.adjust)) +
        geom_point() +
        scale_size_continuous(guide = guide_legend(order = 2))+#, breaks = c( 0.03, 0.07, 0.1, 0.2, 0.3)) +
        scale_color_gradient(low =  "red" , high = "#310A5CFF", name = "adj.P.Value"
                             #,breaks = seq(min(enrich_table$qvalue), max(enrich_table$qvalue), length.out=4),
                             #labels =  round(seq(min(enrich_table$qvalue), max(enrich_table$qvalue), length.out=4),2)
        ) + 
        ylab("Enriched Terms") +
        xlab("-log10(adj.P.Value)") + 
        #theme_bw() + 
        theme_minimal() +
        theme( 
          legend.position="right", 
          #text = element_text(face="bold"),
          axis.text.y = element_text( size = 13), 
          #axis.text = element_text(color = "black", face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) + 
          labs(size = "Hits/Total") +
        theme(axis.text = element_text(color = "black", size = 10, family = "Arial"),
              plot.title = element_text( face = "bold"))+
        ggtitle(paste0("Enrichment Plot of miRNA Targets: ", name))
      
      
      
#print(plot1)
ggsave( paste0(dir,"Enrichment_",  name2, name, ".png"), plot1, dpi = 600, width = 14, height = 12)
}






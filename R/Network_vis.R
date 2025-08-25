setwd(r'{C:\Users\USER\Documents\Github\BC_miRNA}')
dir.create("output")
dir.create("plots")
dir.create("plots/networks/")

#load libraries
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
library(ggraph)
library(tidygraph)


net_df= read.csv("output/mir_targets_de_intersect_subtypes.csv") 
colnames(net_df)= c("Group", "ID", "Target")
head(net_df)
unique(net_df$Group)

for (i in unique(net_df$Group)){
  
      if(i %in% c( "Luminal B", "Luminal A"  )){
        num=1
      }else{
        num=3
      }
  
       name= i
       net= net_df[net_df$Group== i, 2:3 ]
      # get overlapped targets for different miRNAs
      shared_targets <- net %>%
                        group_by(Target) %>%
                        summarise(n_miRNAs = n_distinct(ID)) %>%
                        filter(n_miRNAs >= num)
      
      #subset
      net <- net %>%
             filter(Target %in% shared_targets$Target)
      
      net= net[complete.cases(net) & net$ID != "" & net$Target!= "",]
    
      #for visualization
      g= graph_from_edgelist(as.matrix(net), directed = FALSE)
      # graph for specific nodes
      #graph = induced_graph(graph, vids = V(graph)$name %in% genes)
      # visualize only highest degree nodes
      deg <- igraph::degree(g, mode = "all")
      
      # Step 2: Get top 100 highest-degree nodes
      top_nodes <- names(sort(deg, decreasing = TRUE))[1:100]
      
      # Step 3: Induce subgraph of top 100 nodes
      if (vcount(g) > 100){
      graph <- induced_subgraph(g, vids = top_nodes)
      } else{
      graph= g
      }
      # export
      edge_df <- igraph::as_data_frame(graph, what = "edges")
      write.csv(edge_df, paste0("output/net__", name, ".csv"))
      #vis
      node_colors = rep("grey33", vcount(graph))
      edge_colors = "grey4" # "lavenderblush3" #  #"#F6E8C3"#"grey8"
      node_colors[V(graph)$name %in% net$ID] =  "orange"   #"#180C3CFF" "#F1711FFF"# "#CC4248FF"#"#FFD700", 
      
      
      #ggtable
      graph_tbl = as_tbl_graph(graph) %>%
        activate(nodes) %>%
        mutate(
          # Calculate node degree (frequency of interactions)
          #degree = centrality_degree(mode = "all"),
          # Set size based on degree
          #size = scales::rescale(degree, to = c(2, 10)),
          # col= scales::col_numeric(palette = colors , domain = NULL)(degree),
          # Set color based on degree
          node_colors =node_colors,
          node_shape = "1",
          node_size = 6
        )
      
      
      # Generate the layout and plot
      p= ggraph(graph_tbl, layout = 'circle') + # Try other layouts like 'fr', 'kk' # 'sphere' 'circle'
        geom_edge_link(color = edge_colors, width = 0.08) +
        geom_node_point(aes(color = I(node_colors), size = node_size, shape = node_shape)) +
        geom_node_text(aes(label = name), color = "black", repel = TRUE, size = 3.5) +
        theme_graph() +
        theme(
          plot.title = element_text(hjust = 0.5, color = "black"),
          legend.position = "none"
        ) +
        ggtitle(paste0("Shared Targets of Top-Ranked and DE miRNAs of ",name, " Subtype")) #Regulatory Landscape of Upregulated Driver Genes for TNBC
      
      # Save the plot
      ggsave(paste0("plots/networks/net_targets_intersected_DE_", name, ".png"), p, width = 11, height =12, dpi = 600)

}


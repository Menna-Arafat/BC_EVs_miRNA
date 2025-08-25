setwd(r'{C:\Users\USER\Documents\BC_miRNA_new}')
dir.create("output")
dir.create("plots")

topf= read.csv("sPLSDA_subtypes/output/top_features_component1_subtypes.csv")

for (i in list.files(paste0(getwd(), "/Limma_DE"), pattern = ".csv") ){
    
    name= gsub(".*_([[:alnum:] ]+)\\.csv$", "\\1", i)
    sig= read.csv(paste0("Limma_DE/", i))
    #Intersect top features sPLSDA with DE
    y= intersect(sig$X, topf$X)
    
    write.csv(y,  paste0("output/intersected_sig_" , name, ".csv"), row.names = F)

}

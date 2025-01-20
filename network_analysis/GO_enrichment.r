suppressMessages({
    library(dplyr)
    library(gprofiler2)
})

#import script parameters from config
source("scripts/config/network_analysis/GO_enrichment_config.R")

#import test genes
pset <- read.csv(set_fn, sep='\t', header=TRUE)

#get custom background gene set
pro_ref <- read.csv('data/proteins/protein_annotations_filtered.txt', sep='\t', header=TRUE)
bg <- unique(pro_ref$Protein_name)
    
#run gprofiler2
gostres <- gost(query = pset$Protein_name_A, 
    organism = "hsapiens", 
    custom_bg = bg, 
    correction_method="fdr", 
    evcodes=TRUE, 
    sources=c("GO:MF", "GO:BP", "GO:CC", "KEGG", "REAC", "WP", "MIRNA", "HPA", "HP"))
    

#remove columns, annotate and reorder
go_filt <- subset(gostres$result, select = -c(parents, query, significant))
go_anno <- go_filt[, c("term_id", "term_name", 
        "p_value", "term_size", "query_size", "intersection_size", "precision", "recall", "source", 
        "effective_domain_size", "source_order", "intersection")]

#export results
write.table(go_anno, out_fn, row.names=FALSE, sep="\t", quote = FALSE)

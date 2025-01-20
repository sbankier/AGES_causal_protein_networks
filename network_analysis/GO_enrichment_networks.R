suppressMessages({
    library(dplyr)
    library(gprofiler2)
})

#import script parameters from config
source("scripts/config/network_analysis/GO_enrichment_networks_config.R")

#import protein network
net <- read.csv(net_fn, sep='\t', header=TRUE)

#get custom background gene set
pro_ref <- read.csv('data/proteins/protein_annotations_filtered.txt', sep='\t', header=TRUE)
bg <- unique(pro_ref$Protein_name)

#function to perform GO enrichment using gProfiler
target_go <- function(x) {
    net_tar <- filter(net, S_ID_A == x)
    
    #remove duplicate rows and self edges
    net_un <- net_tar[!duplicated(net_tar[, c("Protein_name_A", "Protein_name_B")]), ]
    net_se <- subset(net_un, net_un$Protein_name_A != net_un$Protein_name_B)
    
    #run gprofiler2
    gostres <- gost(query = net_se$Protein_name_B, 
        organism = "hsapiens", 
        custom_bg = bg, 
        correction_method="fdr", 
        evcodes=TRUE, 
        sources=c("GO:MF", "GO:BP", "GO:CC", "KEGG", "REAC", "WP", "MIRNA", "HPA", "HP"))
    
    return(gostres$result)
}

#get unique A-proteins and perform GO enrichment for all sub-networks
regs <- unique(net$S_ID_A)
go_l <- lapply(regs, target_go)
names(go_l) <- regs
go_df <- bind_rows(go_l, .id = "S_ID_A")

#get the most common terms across all networks
count_terms <- function(df) {
    term_counts <- data.frame(table(df$term_id))
    colnames(term_counts) <- c('term_id', 'number_of_networks')
    df_count <- merge(df, term_counts, by='term_id')
    return(df_count)
}

go_count <- count_terms(go_df)

#remove columns, annotate and reorder
go_count <- subset(go_count, select = -c(parents, query, significant))
go_anno <- merge(go_count, subset(pro_ref, select=c(S_ID, Protein_name)), by.x='S_ID_A', by.y='S_ID')
go_anno <- go_anno[, c("S_ID_A", "Protein_name", "term_id", "term_name", "number_of_networks", 
        "p_value", "term_size", "query_size", "intersection_size", "precision", "recall", "source", 
        "effective_domain_size", "source_order", "intersection")]

#export results
write.table(go_anno, out_fn, row.names=FALSE, sep="\t", quote = FALSE)

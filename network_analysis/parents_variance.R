suppressMessages({
  library(data.table)
  library(dplyr)
})

#import protein expression levels and network
pro <- fread("data/proteins/xpsc_a1_full_adj-sex-age.txt", data.table=FALSE, header=TRUE)
net <- fread('results/findr/filtered/LD_resolution/AGES_findr_networks_adj_1FDR_net10_LD_resolved.tsv', data.table=FALSE)

#import cis-pQTL genotypes and independent SNPs
geno <- fread("data/genotypes/independent/cis_dosage_aproteins_HRC.txt", data.table=FALSE, header=TRUE)
ind <- fread("data/genotypes/independent/independent_cis_signals_aproteins_HRC.csv", data.table=FALSE, header=TRUE)

#format somamer IDs
ind <- ind %>% mutate(somamer = gsub("_", ".", somamer))

#remove all unresolved networks and edges and count number of regulators
net <- net %>% filter(!apply(., 1, function(row) any(grepl(",", row))))
count_b <- net %>% count(S_ID_B, name = "target_frequency")

#assemble as dataframe
net_b <- net %>% 
  inner_join(count_b, by = "S_ID_B") %>% 
  select(S_ID_B, target_frequency) %>% 
  distinct()

#fit all independent SNPs to a linear model to get variance explained for target protein expression
sid_l <- net_b$S_ID_B
r2_l <- numeric(length(sid_l))
r2_adj_l <- numeric(length(sid_l))
for (i in seq_along(sid_l)) {

    #filter protein expression data
    sid <- sid_l[i]
    pro_filt <- select(pro, all_of(c("S_ID", sid)))

    #get somamers of B-protein parents and cis-pQTL (if available)
    par <- filter(net, S_ID_B==sid)$S_ID_A
    par_cis <- c(par, sid)

    #get independent SNPs for all parents
    par_ind <- filter(ind, somamer %in% par_cis)

    #get corresponding genotypes for the independent SNPs
    geno_lab <- c("s_id", par_ind$SNP) %>% unique()
    colnames(geno)[-1] <- sub("_.*", "", colnames(geno)[-1])
    geno_filt <- select(geno, all_of(geno_lab))

    #combine protein and genotypes as single dataframe
    pro_geno <- merge(pro_filt, geno_filt, by.x="S_ID", by.y="s_id") %>% select(-"S_ID")

    #fit to linear model
    formula_str <- paste(colnames(pro_geno)[1], "~", paste(colnames(pro_geno)[-1], collapse = " + "))
    formula_str_bticks <- gsub("([[:alnum:]]+[:.[:alnum:]]+)", "`\\1`", formula_str)
    fit <- lm(formula_str_bticks, data = pro_geno)

    #get variance explained from model
    r2_l[i] <- summary(fit)$r.squared
    r2_adj_l[i] <- summary(fit)$adj.r.squared
}

#format results as dataframe
var_df <- data.frame(
  S_ID = sid_l,
  number_of_regulators = net_b$target_frequency,
  R_squared = r2_l,
  Adjusted_R_squared = r2_adj_l
)

#annotate with protein name
ano <- fread("data/proteins/protein_annotations_filtered.txt", data.table=FALSE, header=TRUE) %>% select("S_ID", "Protein_name")
var_ano <- merge(ano, var_df, by="S_ID")

#export results
fwrite(var_ano, "results/findr/filtered/LD_resolution//AGES_findr_networks_adj_1FDR_net10_LD_targets_varience_explained_by_parents.tsv", sep="\t")

suppressMessages({
  library(data.table)
  library(dplyr)
})

#import protein expression levels and reference
pro <- fread("data/proteins/xpsc_a1_full_adj-sex-age.txt", data.table=FALSE, header=TRUE)

#import cis-pQTL genotypes and independent SNPs
geno <- fread("data/genotypes/independent/cis_dosage_aproteins_HRC.txt", data.table=FALSE, header=TRUE)
ind <- fread("data/genotypes/independent/independent_cis_signals_aproteins_HRC.csv", data.table=FALSE, header=TRUE)

#get unique somamers
sid_l <- unique(ind$somamer)

#fit all independent SNPs to a linear model to get variance explained for protein expression
r2_l <- numeric(length(sid_l))
r2_adj_l <- numeric(length(sid_l))
n_sid_l <- numeric(length(sid_l))
for (i in seq_along(sid_l)) {

  #get independent SNPs for protein
  sid <- sid_l[i]
  ind_filt <- filter(ind, somamer==sid)

  #select protein
  psid <- gsub("_", ".", sid)
  n_sid_l[i] <- psid
  pro_filt <- select(pro, all_of(c("S_ID", psid)))

  #select independent genotypes for corresponding protein
  geno_lab <- c("s_id", ind_filt$SNP)
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
  S_ID = n_sid_l,
  R_squared = r2_l,
  Adjusted_R_squared = r2_adj_l
)

#annotate with protein name
ano <- fread("data/proteins/protein_annotations_filtered.txt", data.table=FALSE, header=TRUE) %>% select("S_ID", "Protein_name")
var_ano <- merge(ano, var_df, by="S_ID")

#export results
fwrite(var_ano, "results/pQTLs/independent_cis_signals_variance_explained_Aproteins_HRC.tsv", sep="\t")

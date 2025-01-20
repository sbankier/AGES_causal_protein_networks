library(data.table)
library(dplyr)
library(coloc)

#set analysis parameters
thrd <- 150000

#import network results
net_mi <- fread("results/colocalisation/AGES_networks_1FDR_net10_MI_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_hf <- fread("results/colocalisation/AGES_networks_1FDR_net10_HF_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_mes <- fread("results/colocalisation/AGES_networks_1FDR_net10_MetS_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_t2d <- fread("results/colocalisation/AGES_networks_1FDR_net10_T2D_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_cac <- fread("results/colocalisation/AGES_networks_1FDR_net10_CAC_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_plq <- fread("results/colocalisation/AGES_networks_1FDR_net10_plaque_overlap.tsv", check.names=TRUE, data.table=FALSE)

#import SNP and protein reference files
pro_ref <- fread("data/proteins/protein_annotations_filtered.txt", check.names=TRUE, data.table=FALSE)
pos_ref <- fread("data/proteins/Soma7000_gene_pos_full.txt", check.names=TRUE, data.table=FALSE)
snp_ref <- fread("data/genotypes/snpinfo_merge.tsv", check.names=TRUE, data.table=FALSE)

#import AGES cis-pQTLs
pq <- fread("results/pQTLs/AGES_cispQTLs_FDR.tsv", check.names=TRUE, data.table=FALSE)

#import GWAS summary statistics
mi <- fread("data/GWAS/MI/GCST011365_buildGRCh37_rsID.tsv", 
    check.names=TRUE, data.table=FALSE, select=c("rsID", "chromosome", "base_pair_location", "beta", "standard_error"))
hf <- fread("data/GWAS/HF/HF_GWAS_rsID.tsv"
    , check.names=TRUE, data.table=FALSE, select=c("rsID", "chr", "pos", "beta", "se"))
mes <- fread("data/GWAS/MetS/Lind2019/UKBB_MetS_alla_Stefan_rsID.tsv"
    , check.names=TRUE, data.table=FALSE, select=c("rsID", "chr", "pos", "beta", "se"))
t2d <- fread("data/GWAS/DIAMANTE/T2D_GWAS_rsID.tsv"
    , check.names=TRUE, data.table=FALSE, select=c("rsID", "chr", "pos", "beta", "se"))
cac <- fread("data/GWAS/CAC/GCST90278456_rsID.tsv"
    , check.names=TRUE, data.table=FALSE, select=c("rsID", "chromosome", "base_pair_location", "beta", "standard_error"))
plq <- fread("data/GWAS/SumStat_plaque/Plaque_meta_032218_rsID.tsv"
    , check.names=TRUE, data.table=FALSE, select=c("rsID", "CHR", "BP", "Effect", "StdErr"))

#set new column names
cn <- c("rsID", "chr", "POS", "beta", "se")
colnames(mi) <- cn
colnames(hf) <- cn
colnames(mes) <- cn
colnames(t2d) <- cn
colnames(cac) <- cn
colnames(plq) <- cn

#filter cis-pQTL summary statistics for network A-proteins
aprot_mi <- unique(net_mi$S_ID)
aprot_hf <- unique(net_hf$S_ID)
aprot_mes <- unique(net_mes$S_ID)
aprot_t2d <- unique(net_t2d$S_ID)
aprot_cac <- unique(net_cac$S_ID)
aprot_plq <- unique(net_plq$S_ID)

run_coloc <- function(sid, gw, n) {

    #get protein reference information
    ref_t<- filter(pro_ref, S_ID == sid)
    pos_t <- filter(pos_ref, SEQ_ID == ref_t[1,2]) %>% slice(1)

    #select SNPs within defined threshold of cis-protein
    up <- pos_t$end+thrd
    dn <- pos_t$start-thrd
    pq_filt <- filter(pq, SEQ_ID == ref_t[1,2] & X.CHROM == pos_t[1,8] & POS < up & POS > dn)
    gw_filt <- filter(gw, chr == pos_t[1,8] & POS < up & POS > dn)

    #calculate beta variance as square of standard error
    pq_filt$betavar <- pq_filt$SE^2
    gw_filt$betavar <- gw_filt$se^2

    #obtain MAF
    maf <- subset(snp_ref, select=c(rsID, MAF))
    pq_maf <- merge(pq_filt, maf, by="rsID")
    gw_maf <- merge(gw_filt, maf, by="rsID")

    #select for coloc input
    pq_sel <- subset(pq_maf, select = c(BETA, betavar, rsID, POS, MAF))
    gw_sel <- subset(gw_maf, select = c(beta, betavar, rsID, POS, MAF))

    #remove duplicate SNPs and missing values
    pq_na <- pq_sel[!duplicated(pq_sel$rsID), ]
    gw_na <- gw_sel[!duplicated(gw_sel$rsID), ]
    pq_na <- na.omit(pq_na)
    gw_na <- na.omit(gw_na)

    #obtain and filter for common SNPs
    com_snps <- intersect(pq_na$rsID, gw_na$rsID)
    pq_stats <- filter(pq_na, rsID %in% com_snps)
    gw_stats <- filter(gw_na, rsID %in% com_snps)

    #obtain data lists for input to coloc
    pq_data <- list(beta=pq_stats$BETA, varbeta=pq_stats$betavar, N=5337, type="quant", snp=pq_stats$rsID, position=pq_stats$POS, MAF=pq_stats$MAF)
    gw_data <- list(beta=gw_stats$beta, varbeta=gw_stats$betavar, N=n, type="quant", snp=gw_stats$rsID, position=gw_stats$POS, MAF=gw_stats$MAF)

    #run coloc for single variant assumption with sensitivity analysis
    coloc_res <- coloc.abf(dataset1=pq_data, dataset2=gw_data)

    return(coloc_res$summary)

}

#Run coloc for all genes in dataset
coloc_mi <- sapply(aprot_mi, run_coloc, gw=mi, n=831000)
coloc_hf <- sapply(aprot_hf, run_coloc, gw=hf, n=977323)
coloc_mes <- sapply(aprot_mes, run_coloc, gw=mes, n=291107)
coloc_t2d <- sapply(aprot_t2d, run_coloc, gw=t2d, n=298957)
coloc_cac <- sapply(aprot_cac, run_coloc, gw=cac, n=26909)
coloc_plq <- sapply(aprot_plq, run_coloc, gw=plq, n=48434)


#format and annotate results
format_res <- function(coloc_res) {
    #assemble as dataframe
    if (class(coloc_res)[1] == "list") {
        coloc_df <- data.frame(do.call(rbind, coloc_res))
    } else {
        coloc_df <- data.frame(t(coloc_res))
    }

    #annotate output
    coloc_df$S_ID <- rownames(coloc_df)
    ano <- subset(pro_ref, select=c(S_ID, Protein_name))
    coloc_ano <- merge(coloc_df, ano, by="S_ID")

    return(coloc_ano)

}

mi_out <- format_res(coloc_mi)
hf_out <- format_res(coloc_hf)
mes_out <- format_res(coloc_mes)
t2d_out <- format_res(coloc_t2d)
cac_out <- format_res(coloc_cac)
plq_out <- format_res(coloc_plq)

#export output
write.table(mi_out, file="results/colocalisation/AGES_networks_1FDR_net10_MI_coloc_res.tsv", quote=FALSE, sep="\t")
write.table(hf_out, file="results/colocalisation/AGES_networks_1FDR_net10_HF_coloc_res.tsv", quote=FALSE, sep="\t")
write.table(mes_out, file="results/colocalisation/AGES_networks_1FDR_net10_MetS_coloc_res.tsv", quote=FALSE, sep="\t")
write.table(t2d_out, file="results/colocalisation/AGES_networks_1FDR_net10_T2D_coloc_res.tsv", quote=FALSE, sep="\t")
write.table(cac_out, file="results/colocalisation/AGES_networks_1FDR_net10_CAC_coloc_res.tsv", quote=FALSE, sep="\t")
write.table(plq_out, file="results/colocalisation/AGES_networks_1FDR_net10_plaque_coloc_res.tsv", quote=FALSE, sep="\t")

library(data.table)
library(dplyr)
library(coloc)
library(susieR)

#import network results
net_mi <- fread("results/colocalisation/AGES_networks_1FDR_net10_MI_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_hf <- fread("results/colocalisation/AGES_networks_1FDR_net10_HF_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_mes <- fread("results/colocalisation/AGES_networks_1FDR_net10_MetS_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_t2d <- fread("results/colocalisation/AGES_networks_1FDR_net10_T2D_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_cac <- fread("results/colocalisation/AGES_networks_1FDR_net10_CAC_overlap.tsv", check.names=TRUE, data.table=FALSE)
net_plq <- fread("results/colocalisation/AGES_networks_1FDR_net10_plaque_overlap.tsv", check.names=TRUE, data.table=FALSE)

#Get network proteins with LD matrices
ld_net <- fread("results/colocalisation/AGES_findr_networks_adj_1FDR_net10_aprot.tsv", check.names=TRUE, data.table=FALSE)

#import SNP and protein reference files
pro_ref <- fread("data/proteins/protein_annotations_filtered.txt", check.names=TRUE, data.table=FALSE)
pos_ref <- fread("data/proteins/Soma7000_gene_pos_full.txt", check.names=TRUE, data.table=FALSE)
snp_ref <- fread("data/genotypes/snpinfo_merge.tsv", check.names=TRUE, data.table=FALSE)
ld_ref <- fread("V:/LD_matrices/pos_id_mapping.csv", check.names=TRUE, data.table=FALSE)

#import AGES cis-pQTL results
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

#function to run susie coloc for overlap between pQTL and GWAS trait
run_coloc_susie <- function(label, net, gw, thrd, n, L, max_iter, gwas_harm) {

    #set new column names for GWAS results
    cn <- c("rsID", "chr", "POS", "beta", "se")
    colnames(gw) <- cn

    #filter cis-pQTL summary statistics for network A-proteins
    ld_net_gw <- merge(ld_net, net, by.x="S_ID_A", by.y="S_ID")
    aprot <- unique(ld_net_gw$S_ID)

    #run for each overlapping pQTL
    coloc_res <- list()
    for (sid in aprot) {

        message("Processing locus ", sid, "...")

        #import LD matrix
        ld_fn <- paste0("V:/LD_matrices/", sid, "_LDmatrix.rds")
        ld <- readRDS(ld_fn)

        #protein reference info
        ref_t <- filter(pro_ref, S_ID == sid)
        pos_t <- filter(pos_ref, SEQ_ID == ref_t[1,2]) %>% slice(1)

        #SNPs within cis-window
        up <- pos_t$end + thrd
        dn <- pos_t$start - thrd
        pq_filt <- filter(pq, SEQ_ID == ref_t[1,2] & X.CHROM == pos_t[1,8] & POS < up & POS > dn)
        gw_filt <- filter(gw, chr == pos_t[1,8] & POS < up & POS > dn)

        #beta variance
        pq_filt$betavar <- pq_filt$SE^2
        gw_filt$betavar <- gw_filt$se^2

        #merge with MAF
        maf <- subset(snp_ref, select=c(rsID, MAF))
        pq_maf <- merge(pq_filt, maf, by="rsID")
        gw_maf <- merge(gw_filt, maf, by="rsID")

        #annotate with GRCh38 markers
        ld_ano <- subset(ld_ref, select=c(rsid, id))
        pq_ld <- merge(pq_maf, ld_ano, by.x="rsID", by.y="rsid")
        gw_ld <- merge(gw_maf, ld_ano, by.x="rsID", by.y="rsid")

        #select for coloc input
        pq_sel <- subset(pq_ld, select = c(BETA, betavar, id, POS, MAF))
        gw_sel <- subset(gw_ld, select = c(beta, betavar, id, POS, MAF))

        #remove duplicates and NAs
        pq_na <- na.omit(pq_sel[!duplicated(pq_sel$id), ])
        gw_na <- na.omit(gw_sel[!duplicated(gw_sel$id), ])

        #common SNPs
        com_snps <- intersect(pq_na$id, gw_na$id)
        pq_stats <- filter(pq_na, id %in% com_snps)
        gw_stats <- filter(gw_na, id %in% com_snps)

        #harmonise GWAS effect direction if required
        if (isTRUE(gwas_harm)) {
        gw_stats$beta <- -gw_stats$beta
        }

        #filter LD matrix for common SNPs
        test_snps <- pq_stats$id
        ld_filt <- ld[test_snps, test_snps]

        #prepare data lists for coloc
        pq_data <- list(beta=pq_stats$BETA, varbeta=pq_stats$betavar, N=5337,
                        type="quant", snp=pq_stats$id, position=pq_stats$POS, MAF=pq_stats$MAF, LD=ld_filt)
        gw_data <- list(beta=gw_stats$beta, varbeta=gw_stats$betavar, N=n,
                        type="quant", snp=gw_stats$id, position=gw_stats$POS, MAF=gw_stats$MAF, LD=ld_filt)

        #fit SuSiE for pQTL
        susie_pq <- try(
            susie_rss(
            bhat=pq_data$beta,
            shat=sqrt(pq_data$varbeta),
            R=ld_filt,
            L=L,
            n=pq_data$N,
            max_iter=max_iter,
            estimate_residual_variance=FALSE
            ), silent=TRUE
        )
        if (inherits(susie_pq, "try-error") || !isTRUE(susie_pq$converged)) {
            message("Skipping locus ", sid, " pQTL SuSiE did not converge")
            next
        }

        #fit SuSiE for GWAS
        susie_gw <- try(
            susie_rss(
            bhat=gw_data$beta,
            shat=sqrt(gw_data$varbeta),
            R=ld_filt,
            L=L,
            n=gw_data$N,
            max_iter = max_iter,
            estimate_residual_variance=FALSE
            ), silent=TRUE
        )
        if (inherits(susie_gw, "try-error") || !isTRUE(susie_gw$converged)) {
            message("Skipping locus ", sid, " GWAS SuSiE did not converge")
            next
        }

        #run coloc if susie results available for both pQTL and GWAS
        if (requireNamespace("susieR", quietly = TRUE)) {
            susie.res <- coloc.susie(susie_pq, susie_gw)
            susie_out <- susie.res$summary
            susie_out$S_ID <- sid
            message(sid, " complete.")
            coloc_res[[sid]] <- susie_out
        }
    }

    #combine results into a single dataframe
    coloc_res_df <- rbindlist(
    coloc_res[!sapply(coloc_res, is.null)], 
    use.names = TRUE, 
    fill = TRUE)

    #get protein annotations
    pro_ref <- ld_net_gw %>% 
    select(S_ID = S_ID_A, Protein_name = Protein_name_A) %>% 
    distinct()

    #annotate coloc results
    coloc_res_ano <- coloc_res_df %>%
    left_join(pro_ref, by = "S_ID", relationship = "many-to-many") %>%
    select(S_ID, Protein_name, everything())
    coloc_res_ano$trait <- label

    #export results
    out_fn <- paste0("results/colocalisation/susie/AGES_networks_1FDR_net10_", label, "_coloc_susie_res.tsv")
    fwrite(coloc_res_ano, out_fn, sep="\t")

    return(coloc_res_ano)
}

#run coloc susie for all GWAS traits
mi_coloc <- run_coloc_susie(label="MI", net=net_mi, gw=mi, thrd=150000, n=831000, L=5, max_iter=5000, gwas_harm=FALSE)
hf_coloc <- run_coloc_susie(label="HF", net=net_hf, gw=hf, thrd=150000, n=977323, L=5, max_iter=5000, gwas_harm=TRUE)
mes_coloc <- run_coloc_susie(label="MES", net=net_mes, gw=mes, thrd=150000, n=291107, L=5, max_iter=5000, gwas_harm=FALSE)
t2d_coloc <- run_coloc_susie(label="T2D", net=net_t2d, gw=t2d, thrd=150000, n=298957, L=5, max_iter=5000, gwas_harm=TRUE)
cac_coloc <- run_coloc_susie(label="CAC", net=net_cac, gw=cac, thrd=150000, n=26909, L=5, max_iter=5000, gwas_harm=FALSE)
plq_coloc <- run_coloc_susie(label="plaque", net=net_plq, gw=plq, thrd=150000, n=48434, L=5, max_iter=5000, gwas_harm=TRUE)

#combine as single dataframe
coloc_res_full <- list(mi_coloc, hf_coloc, mes_coloc, t2d_coloc, cac_coloc, plq_coloc)
coloc_res_full_df <- bind_rows(coloc_res_full)

#export results
coloc_full_fn <- "results/colocalisation/susie/AGES_networks_1FDR_net10_full_coloc_susie_res.tsv"
fwrite(coloc_res_full_df, coloc_full_fn, sep="\t")

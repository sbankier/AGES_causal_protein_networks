import pandas as pd

#import network results
net=pd.read_csv('results/findr/filtered/AGES_findr_networks_adj_1FDR_net10.tsv', delimiter='\t')
aprot=net[['S_ID_A']].drop_duplicates()
aprot.columns=['S_ID']

#import AGES cis-pQTLs and SNP reference
pq=pd.read_csv('results/pQTLs/AGES_cispQTLs_FDR.tsv', delimiter='\t')
snp_ref=pd.read_csv('data/genotypes/snpinfo_merge.tsv', delimiter='\t')

#import GWAS summary statistics
mi=pd.read_csv('data/GWAS/MI/GCST011365_buildGRCh37.tsv', delimiter='\t')
hf=pd.read_csv('data/GWAS/HF/HF_GWAS.txt', delimiter='\t')
mes=pd.read_csv('data/GWAS/MetS/Lind2019/UKBB_MetS_alla_Stefan.txt', delimiter='\t')
t2d=pd.read_csv('data/GWAS/DIAMANTE/T2D_GWAS.txt', delimiter='\t')
cac=pd.read_csv('data/GWAS/CAC/GCST90278456.tsv', delimiter='\t')
plq=pd.read_csv('data/GWAS/SumStat_plaque/Plaque_meta_032218.csv', delimiter=',')

#add rsIDs to GWAS results
mi_rs=pd.merge(mi, snp_ref[['rsID', 'chr', 'pos']], left_on=['chromosome', 'base_pair_location'], right_on=['chr', 'pos'], how='inner').drop('variant_id', axis=1)
hf_rs=pd.merge(hf, snp_ref[['rsID', 'chr', 'pos']], on=['chr', 'pos'], how='inner')
mes_rs=pd.merge(mes, snp_ref[['rsID', 'chr', 'pos']], on=['chr', 'pos'], how='inner').drop(['rsid'], axis=1)
t2d_rs=pd.merge(t2d, snp_ref[['rsID', 'chr', 'pos']], on=['chr', 'pos'], how='inner')
cac_rs=pd.merge(cac, snp_ref[['rsID', 'chr', 'pos']], left_on=['chromosome', 'base_pair_location'], right_on=['chr', 'pos'], how='inner').drop(['variant_id'], axis=1)
plq_rs=pd.merge(plq, snp_ref[['rsID', 'chr', 'pos']], left_on=['CHR', 'BP'], right_on=['chr', 'pos'], how='inner').drop(['ID'], axis=1)

#export rsID annotated GWAS results
mi_rs.to_csv('data/GWAS/MI/GCST011365_buildGRCh37_rsID.tsv', sep='\t', index=None)
hf_rs.to_csv('data/GWAS/HF/HF_GWAS_rsID.tsv', sep='\t', index=None)
mes_rs.to_csv('data/GWAS/MetS/Lind2019/UKBB_MetS_alla_Stefan_rsID.tsv', sep='\t', index=None)
t2d_rs.to_csv('data/GWAS/DIAMANTE/T2D_GWAS_rsID.tsv', sep='\t', index=None)
cac_rs.to_csv('data/GWAS/CAC/GCST90278456_rsID.tsv', sep='\t', index=None)
plq_rs.to_csv('data/GWAS/SumStat_plaque/Plaque_meta_032218_rsID.tsv', sep='\t', index=None)

#format GWAS results
mi_filt=mi_rs[['rsID', 'beta', 'standard_error', 'p_value']]
hf_filt=hf_rs[['rsID', 'beta', 'se', 'pval']]
mes_filt=mes_rs[['rsID', 'beta', 'se', 'p_value']]
t2d_filt=t2d_rs[['rsID', 'beta', 'se', 'pval']]
cac_filt=cac_rs[['rsID', 'beta', 'standard_error', 'p_value']]
plq_filt=plq_rs[['rsID', 'Effect', 'StdErr', 'P_value']]

#filter combined results for A-proteins
def com_res(gw):
    gw.columns=['rsID', 'beta', 'se', 'pv']
    pq_a=pd.merge(pq, aprot, on='S_ID', how='inner')
    pq_gw=pd.merge(pq_a, gw, on='rsID', how='inner')
    pq_gw_sig=pq_gw.loc[(pq_gw['q-value'] <= 0.05) & (pq_gw['pv'] <= 5e-5)]
    
    return pq_gw_sig

#obtain combined pQTL/ GWAS results
pq_mi_sig=com_res(mi_filt)
pq_hf_sig=com_res(hf_filt)
pq_mes_sig=com_res(mes_filt)
pq_t2d_sig=com_res(t2d_filt)
pq_cac_sig=com_res(cac_filt)
pq_plq_sig=com_res(plq_filt)

#export combined filtered results
pq_mi_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_MI_overlap.tsv', sep='\t', index=None)
pq_hf_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_HF_overlap.tsv', sep='\t', index=None)
pq_mes_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_MetS_overlap.tsv', sep='\t', index=None)
pq_t2d_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_T2D_overlap.tsv', sep='\t', index=None)
pq_cac_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_CAC_overlap.tsv', sep='\t', index=None)
pq_plq_sig.to_csv('results/colocalisation/AGES_networks_1FDR_net10_plaque_overlap.tsv', sep='\t', index=None)
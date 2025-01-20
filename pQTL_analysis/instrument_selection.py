import pandas as pd
from multipy import fdr

#import cis-pQTLs
cis=pd.read_csv('results/pQTLs/cispQTLs.txt', delimiter='\t')

#get p-values from -log10pv
pvl=[]
for p in cis['LOG10_P']:
	pv=pow(10, -p)
	pvl.append(pv)
cis['p-value']=pvl

#obtain FDR adjusted q-values
qv=fdr.qvalue(cis['p-value'].values)
cis['q-value']=qv[1]

#add rsID
snp_ref=pd.read_csv('data/genotypes/snpinfo_merge.tsv', delimiter='\t')
cis_ref=pd.merge(cis, snp_ref[['SNP', 'rsID']], left_on='ID', right_on='SNP', how='inner').drop('SNP', axis=1)

#annotate SOMA IDs
pro_ref=pd.read_csv('data/proteins/protein_annotations_filtered.txt', delimiter='\t')
cis_ref=pd.merge(cis_ref, pro_ref[['SeqId', 'S_ID', 'Target']], left_on='SEQ_ID', right_on='SeqId', how='inner').drop('SeqId', axis=1)

#reorientate columsn and export
cis_ref=cis_ref[['rsID', 'SEQ_ID', 'S_ID', 'Target', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 
		'A1', 'TEST', 'OBS_CT', 'BETA','SE', 'T_STAT', 'LOG10_P', 'p-value', 'q-value']]

cis_ref.to_csv('results/pQTLs/AGES_cispQTLs_FDR.tsv', sep='\t', index=None)

#filter for top pQTLs given a 5% FDR threshold 
cis_fdr=cis_ref.loc[cis_ref['q-value'] <= 0.05]
cis_ins=cis_fdr.sort_values('LOG10_P', ascending=False).drop_duplicates('SEQ_ID', keep='first')

#export pQTL instrumnets
cis_ins.to_csv('results/pQTLs/AGES_cispQTLs_instrument_5FDR.tsv', sep='\t', index=None)
import pandas as pd
import numpy as np
import yaml
import math
import findr

"""Iterates aross list of proteins with best cis-pQTL and calculates all 5 tests from Findr 
using best eQTLs as causal anchors. Returns output of all Findr tests and manual calculation 
of CIT and P2*P5 test combination"""

#read config
config_path='scripts/config/findr_cis-anchor.yml'

with open(config_path, 'r') as config_file:
        config = yaml.safe_load(config_file)

#import protein expression and genotype data and protein reference
pro_fn=config['paths']['pro_fn']
pro=pd.read_csv(pro_fn, delimiter='\t')
geno=pd.read_csv('data/genotypes/AGES_cispQTLs_instrument_5FDR_dosages_rsID_full.txt', delimiter='\t')
pro_ref=pd.read_csv('data/proteins/protein_annotations_filtered.txt', delimiter='\t')

#import cis-pQTL instruments for causal inference
pqtl=pd.read_csv('results/pQTLs/AGES_cispQTLs_instrument_5FDR_full.tsv', delimiter='\t')

#obtain common samples between genotype and proteins
com=pd.merge(pro[['S_ID']], geno[['rsID']], left_on='S_ID', right_on='rsID', how='inner')
pro_com=pd.merge(com[['S_ID']], pro, on='S_ID', how='inner').drop('PROJ', axis=1)
geno_com=pd.merge(com[['rsID']], geno, on='rsID', how='inner')

#reorientate dataframes
def flip_df(df):
	df=df.T.reset_index()
	df.columns = df.iloc[0]
	df=df.iloc[1:]
	return df

prof=flip_df(pro_com)
genof=flip_df(geno_com)

#get instrument genotypes and proteins
pq_geno_all=pd.merge(pqtl[['rsID', 'S_ID']], genof, on='rsID', how='inner')
pq_pro_all=pd.merge(pqtl[['S_ID']], prof, on='S_ID', how='inner')

#iterate across instruments and run findr each time for A-proteins against all B proteins
findr_l=[]
for pq in pqtl['S_ID']:

	#get A-protein expression and genotype
	pq_geno=pq_geno_all.loc[pq_geno_all['S_ID'] == pq]
	pq_pro=pq_pro_all.loc[pq_pro_all['S_ID'] == pq]

	#move A protein to top of B matrix for nodiag
	pro_filt=prof.loc[prof['S_ID'] != pq]
	pro_sort=pd.concat([pq_pro, pro_filt])

	#remove columns and convert to format for findr
	dg=pq_geno.drop(['rsID', 'S_ID'], axis=1).values.astype(np.uint8)
	d=pq_pro.drop('S_ID', axis=1).values.astype(np.float32)
	dt=pro_sort.drop(['S_ID'], axis=1).values.astype(np.float32)

	#raise error if input shapes are different
	if np.shape(dg) != np.shape(d):
		with open('results/findr_error.log', 'a') as e_out:
			e_out.write(pq+' failed: skipping...'+'\n')
			e_out.close()
		print(pq+' failed: skipping...')

		continue
	else:
		pass

	#run findr and obtain output of all tests
	l=findr.lib(path="/home/sean/temp/findr/libfindr.so",loglv=6,rs=0,nth=0)
	ans=l.pijs_gassist(dg,d,dt,na=None,nodiag=True)

	#convert from p matrices to long
	def wide_t_long(p):
		df=pd.DataFrame(ans[p], index=pq_pro['S_ID'].values, columns=pro_sort['S_ID'].values)
		plong=df.stack().reset_index()
		plong.columns=['S_ID_A', 'S_ID_B', p]
		return plong[['S_ID_A', 'S_ID_B']], plong[[p]]

	#annotate S_IDs with protein target
	findr_df=wide_t_long('p2')[0]
	findr_df=pd.merge(findr_df, pro_ref[['S_ID', 'Target']], left_on='S_ID_A', right_on='S_ID', how='left').drop('S_ID', axis=1)
	findr_df=pd.merge(findr_df, pro_ref[['S_ID', 'Target']], left_on='S_ID_B', right_on='S_ID', how='left').drop('S_ID', axis=1)
	findr_df.columns=['S_ID_A', 'S_ID_B', 'Target_A', 'Target_B']
	findr_df=findr_df[['S_ID_A', 'Target_A', 'S_ID_B', 'Target_B']]

	#obtain dataframe with all tests
	for p in ['p2', 'p3', 'p4', 'p5']:
		ldf=wide_t_long(p)[1]
		findr_df[p]=ldf

	#calculate test combinations
	p2p3=findr_df['p2']*findr_df['p3']
	p2p5=findr_df['p2']*findr_df['p5']
	findr_df['p2p3']=p2p3
	findr_df['p2p5']=p2p5

	#create list of dataframes for all results
	findr_l.append(findr_df)

#export results
out=pd.concat(findr_l)
out_fn=config['paths']['out_fn']
out.to_csv(out_fn, sep='\t', index=None)
import pandas as pd

#import genotypes and SNP reference
geno=pd.read_csv('data/genotypes/AGES_cispQTLs_instrument_5FDR_dosages_full.txt', delimiter='\t')
ref=pd.read_csv('data/genotypes/snpinfo_merge.tsv', delimiter='\t')

#flip dataframe
def flip_df(df):
	df=df.T.reset_index()
	df.columns = df.iloc[0]
	df=df.iloc[1:]
	return df

genof=flip_df(geno)

#re-annotate SNP ID
idl=[]
for x in genof['IID']:
	s=x.split('_')
	idl.append(s[0])

genof['IID']=idl

#annotate genotypes with rsIDs
geno_rs=pd.merge(ref[['SNP', 'rsID']], genof, left_on='SNP', right_on='IID', 
		how='inner').drop(['SNP', 'IID'], axis=1)

#reorientate and export
geno_out=flip_df(geno_rs)
geno_out.to_csv('data/genotypes/AGES_cispQTLs_instrument_5FDR_dosages_rsID_full.txt', sep='\t', index=None)
import pandas as pd
pd.options.mode.chained_assignment = None
from functools import reduce
from multipy import fdr

#import trait reference
t_ref=pd.read_csv('results/pheno_assoc/trait_reference.txt', delimiter='\t')
mi_ref=t_ref.loc[t_ref['MI_trait'] == True]

#import network eigneprotein phenotype associations and filter for 5% FDR
eign_ph_l=mi_ref.loc[mi_ref['eigenprotein'] == True]['label'].tolist()
eign_dfl=[]
for lab in eign_ph_l:
    df_eign=pd.read_csv('results/pheno_assoc/eigenproteins/PC1_associations_'+lab+'.txt', delimiter='\t')
    lab_n=lab+'_5FDR'
    df_eign[lab_n] = (df_eign['P-value'] < 0.00027).astype(int)
    eign_dfl.append(df_eign[['Network', lab_n]])

#export boolean dataframe for eigenprotein associations
eign_sec=reduce(lambda df1,df2: pd.merge(df1,df2,on='Network', how='inner'), eign_dfl)
eign_sec.to_csv('results/pheno_assoc/eigenproteins/PC1_associations_5FDR_intersection.tsv', sep='\t', index=None)

#import network A-protein phenotype associations and filter for 5% FDR
aprot_ph_l=mi_ref.loc[mi_ref['A-protein'] == True]['label'].tolist()
aprot_dfl=[]
for lab in aprot_ph_l:
    df_aprot=pd.read_csv('results/pheno_assoc/A-proteins/A_protein_associations_'+lab+'.txt', delimiter='\t')
    lab_n=lab+'_5FDR'
    df_aprot[lab_n] = (df_aprot['P-value adj'] < 0.05).astype(int)
    aprot_dfl.append(df_aprot[['A-Protein', lab_n]])

#export boolean dataframe for A-protein associations
aprot_sec = reduce(lambda df1,df2: pd.merge(df1,df2,on='A-Protein', how='inner'), aprot_dfl)
aprot_sec.to_csv('results/pheno_assoc/A-proteins/A_protein_associations_5FDR_intersection.tsv', sep='\t', index=None)
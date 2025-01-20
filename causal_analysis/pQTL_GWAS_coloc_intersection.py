import pandas as pd
pd.options.mode.chained_assignment = None
from functools import reduce

#phenotype values
ph=['MI', 'HF', 'MetS', 'T2D', 'CAC', 'Plaque']

#combine results in single dataframe
coloc_l=[]
for x in ph:
    coloc=pd.read_csv('results/colocalisation/AGES_networks_1FDR_net10_'+x+'_coloc_res.tsv', delimiter='\t')
    coloc['trait']=x
    coloc_l.append(coloc)

#export combined results
coloc_all=pd.concat(coloc_l)
coloc_all.to_csv('results/colocalisation/AGES_networks_1FDR_net10_total_coloc_res.tsv', sep='\t', index=None)

#Set PP threshold for coloc results and find intersection between phenotype associations
def coloc_filt(pp):
    coloc_res=[]
    for x in ph:
        coloc=pd.read_csv('results/colocalisation/AGES_networks_1FDR_net10_'+x+'_coloc_res.tsv', delimiter='\t')
        coloc_un=coloc.sort_values('PP.H4.abf', ascending=False).drop_duplicates('Protein_name')
        coloc_un[x]=(coloc_un['PP.H4.abf'] >= pp).astype(int)
        coloc_res.append(coloc_un[['S_ID', 'Protein_name', x]])
    
    #combine as single dataframe
    coloc_sig=reduce(lambda df1,df2: pd.merge(df1,df2,on=['S_ID', 'Protein_name'], how='outer'), coloc_res).fillna(0)
    
    #compress values
    coloc_com=coloc_sig.drop(columns=[col for col in coloc_sig.columns if coloc_sig[col].eq(0).all()])
    coloc_com=coloc_com.groupby('Protein_name', as_index=False).max()

    #obtain frequency of colocalisations
    nph=coloc_com.drop(['Protein_name', 'S_ID'], axis=1).columns
    freq=coloc_com[nph].sum(axis=1)
    coloc_com['frequency']=freq
    coloc_out=coloc_com.loc[coloc_com['frequency'] >= 1].drop('S_ID', axis=1).sort_values('frequency', ascending=False)

    return(coloc_out)

#obtain colocalisations at different thesholds
coloc_5=coloc_filt(0.5)
coloc_9=coloc_filt(0.9)

#export results
coloc_5.to_csv('results/colocalisation/AGES_networks_1FDR_net10_coloc_res_intersection_0.5.tsv', sep='\t', index=None)
coloc_9.to_csv('results/colocalisation/AGES_networks_1FDR_net10_coloc_res_intersection_0.9.tsv', sep='\t', index=None)
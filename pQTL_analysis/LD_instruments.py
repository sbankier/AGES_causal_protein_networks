import pandas as pd

#import cis-pQTL instrumnets
pq=pd.read_csv('results/pQTLs/AGES_cispQTLs_instrument_5FDR_full.tsv', delimiter='\t')

#import reconstructed network
net=pd.read_csv('results/findr/filtered/AGES_findr_networks_adj_1FDR_net10.tsv', delimiter='\t')
net_un=net.sort_values('number_of_targets', ascending=False)[['S_ID_A', 'Protein_name_A']].drop_duplicates('Protein_name_A')
net_un.columns=['S_ID', 'Protein_name']

#filter for network proteins
pq_net=pd.merge(pq, net_un, on='S_ID', how='inner')
pq_un=pq_net['rsID'].unique().tolist()

#import genotypes
geno=pd.read_csv('data/genotypes/AGES_cispQTLs_instrument_5FDR_dosages_rsID_full.txt', delimiter='\t')

#obtain LD matrix from instruments
pq_geno=geno[pq_un]
ld_mat=pq_geno.corr(method='pearson').reset_index().rename(columns={'index': 'rsID'})

#export LD matrix
ld_mat.to_csv('results/pQTLs/Findr_1FDR_net10_instrument_LD_matrix.tsv', sep='\t', index=None)

#get R2 from LD matrix
ld_mat.iloc[:, 1:]=ld_mat.iloc[:, 1:] ** 2

#obtain blocks of protein with either shared instrument or in high LD
ld_d={}
for rs in pq_un:
    pq_cor=ld_mat[['rsID', rs]]
    pq_cor_sig=pq_cor.loc[pq_cor[rs] >= 0.5][['rsID']]
    pq_ld=pd.merge(pq_net[['rsID', 'S_ID', 'Protein_name']], pq_cor_sig, on='rsID', how='inner')

    #get LD block label
    ld_blk=sorted(pq_ld['Protein_name'].tolist())
    ld_lab=','.join(ld_blk)

    #assemble LD blocks
    if len(pq_ld) > 1:
        ld_d[ld_lab]=pq_ld
    else:
        pass

#combine all LD blocks as single dataframe
ld_net=pd.concat(ld_d).reset_index().drop('level_1', axis=1)
ld_net.rename(columns={'level_0': 'LD_block'}, inplace=True)

#count number of proteins in LD block
ld_count=ld_net['LD_block'].value_counts().to_dict()
ld_countdf=pd.DataFrame(ld_count.items(), columns=['LD_block', 'LD_count'])
ld_out=pd.merge(ld_countdf, ld_net, on='LD_block', how='inner')

#export LD blocks
ld_out.to_csv('results/pQTLs/AGES_cispQTLs_instrument_5FDR_shared_LD.tsv', sep='\t', index=None)
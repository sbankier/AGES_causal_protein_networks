import pandas as pd
from itertools import combinations
from scipy import stats
pd.options.mode.chained_assignment = None

#import reconstructed network
net=pd.read_csv('results/findr/filtered/AGES_findr_networks_adj_1FDR_net10.tsv', delimiter='\t')

#remove duplicate somamers
sid_un=net.sort_values('number_of_targets', ascending=False).drop_duplicates('Protein_name_A')[['S_ID_A']]
net_un=pd.merge(net, sid_un, on='S_ID_A', how='inner').drop(['Target_A', 'Target_B'], axis=1)

#import protein expression data
pro=pd.read_csv('data/proteins/xpsc_a1_full_adj-sex-age.txt', delimiter='\t')

#import LD blocks
ld_bl=pd.read_csv('results/pQTLs/AGES_cispQTLs_instrument_5FDR_shared_LD.tsv', delimiter='\t')

#obtain correlations between proteins within LD blocks
def pro_cor(ld):
    
    #get proteins within LD block
    sh_pq=ld_bl.loc[ld_bl['LD_block'] == ld]['S_ID'].tolist()

    #get all pairwise combinations
    pq_com=list(combinations(sh_pq, 2))
    pq_com_ref=[(ld_bl[ld_bl['S_ID'] == id1]['Protein_name'].values[0], ld_bl[ld_bl['S_ID'] == id2]['Protein_name'].values[0]) 
                for id1, id2 in pq_com]

    #calculate protein-protein correlations
    cor_l=[]
    for s_pair, p_pair in zip(pq_com, pq_com_ref):
        s1, s2 = s_pair
        p1, p2 = p_pair
        sp=stats.spearmanr(pro[s1], pro[s2])
        cor=sp[0]
        pv=sp[1]

        #get intersection/ union as percentage
        s1_net=net_un.loc[net_un['S_ID_A'] == s1].drop_duplicates('Protein_name_B')
        s2_net=net_un.loc[net_un['S_ID_A'] == s2].drop_duplicates('Protein_name_B')
        pair_union=list(set(s1_net['Protein_name_B']).union(set(s2_net['Protein_name_B'])))
        pair_intsc=list(set(s1_net['Protein_name_B']).intersection(set(s2_net['Protein_name_B'])))
        su_per=(len(pair_intsc)/len(pair_union))*100

        #assemble as Dataframe
        cor_df=pd.DataFrame({'S_ID_1': [s1], 'Protein_1': [p1], 'S_ID_2': [s2], 'Protein_2': [p2], 
                            'Pearson_r': [cor],
                            'Pearson_pv': [pv],
                            'Protein_1_network_size': [len(s1_net)],
                            'Protein_2_network_size': [len(s2_net)],
                            'Protein_pair_intersection': [len(pair_intsc)], 
                            'Protein_pair_union': [len(pair_union)],
                            'intersection_union_percentage': [su_per]})
        cor_l.append(cor_df)

    #concatenate and add LD block label
    cor_df_full=pd.concat(cor_l)
    cor_df_full['LD_block']=ld

    return cor_df_full

#get unqiue LD blocks
ld_un=ld_bl['LD_block'].unique()

#obtain single dataframe of pairwise correlations for all LD block proteins
df_l=[]
for ld in ld_un:
    df=pro_cor(ld)
    df_l.append(df)

ld_cor_df=pd.concat(df_l)

#export LD block protein correlations/ target overlap
ld_cor_df.to_csv('results/findr/filtered/LD_resolution/AGES_findr_networks_adj_1FDR_net10_LD_protein_correlations.tsv', sep='\t', index=None)

#flip and concatenate protein pairwise correlations
ld_cor_f=ld_cor_df.copy()
f_cols={'S_ID_1': 'S_ID_2', 'S_ID_2': 'S_ID_1', 'Protein_1': 'Protein_2', 'Protein_2': 'Protein_1', 
        'Protein_1_network_size': 'Protein_2_network_size', 'Protein_2_network_size': 'Protein_1_network_size'}
ld_cor_f.rename(columns=f_cols, inplace=True)
ld_cor_mir=pd.concat([ld_cor_df, ld_cor_f])

#get LD block proteins as single column
ld_pro=pd.DataFrame({'Protein_name': ld_cor_df['Protein_1'].tolist()+ld_cor_df['Protein_2'].tolist(),
                    'S_ID': ld_cor_df['S_ID_1'].tolist()+ld_cor_df['S_ID_2'].tolist(),
                    'LD_block': ld_cor_df['LD_block'].tolist()+ld_cor_df['LD_block'].tolist()}).drop_duplicates().sort_values('LD_block')

#get a summary for each protein within an LD block
pro_n=[]
ld_union=[]
for p in ld_pro['Protein_name']:

    #get network for A-protein
    net_pro=net_un.loc[net_un['Protein_name_A'] == p].drop_duplicates('Protein_name_B')
    pro_n.append(len(net_pro))

    #get union of targets for all A-proteins within LD block
    ld=ld_pro.loc[ld_pro['Protein_name'] == p]['LD_block'].values[0]
    lds=ld.split(',')
    ld_net=net_un.loc[net_un['Protein_name_A'].isin(lds)].drop_duplicates('Protein_name_B')
    ld_union.append(len(ld_net))

ld_pro['A-protein_network_size']=pro_n
ld_pro['LD_network_union_size']=ld_union

#resolve LD blocks as either merged or independent networks
net_res=net_un.copy()
ld_res=[]
for ld in ld_un:

    #identify if LD block contains mutual network targets
    lds=ld.split(',')
    ld_net=net_un.loc[(net_un['Protein_name_A'].isin(lds)) & (net_un['Protein_name_B'].isin(lds))]
    mut_tar=set(ld_net['Protein_name_A'].tolist()+ld_net['Protein_name_B'].tolist())

    #resolve LD blocks where there are mutual targets
    if len(mut_tar) > 0:

        #get union/ intersection percentages for mutal target proteins
        mut_per=ld_cor_mir.loc[ld_cor_mir['Protein_1'].isin(list(mut_tar))]

        #Identify mutual targets where intersection/union percentage is > 60%
        mut_sig=mut_per.loc[mut_per['intersection_union_percentage'] > 60]
        mut_sig_lab=list(set(mut_sig['Protein_1'].tolist()+mut_sig['Protein_2'].tolist()))
        mut_bl=','.join(mut_sig_lab)

        #get SID labels for merged blocks
        mut_sig_lab_sid=list(set(mut_sig['S_ID_1'].tolist()+mut_sig['S_ID_2'].tolist()))
        mut_bl_sid=','.join(mut_sig_lab_sid)

        #replace protein names with updated label for merged networks
        if len(mut_bl) > 0:
            net_res=net_res.replace(to_replace=mut_sig_lab, value=mut_bl)
            net_res=net_res.replace(to_replace=mut_sig_lab_sid, value=mut_bl_sid)
        
        #extract protein names for independent networks
        ind_tar=set(lds)-set(mut_sig_lab)

        #record resolved LD block labels
        if len(ind_tar) > 0 and len(mut_bl) > 0:
            ld_ind_lab=' + '.join(list(ind_tar))
            ld_res.append(mut_bl+' + '+ld_ind_lab)
        elif len(ind_tar) > 0 and len(mut_bl) == 0:
            ld_ind_lab=' + '.join(list(ind_tar))
            ld_res.append(ld_ind_lab)
        elif len(ind_tar) == 0 and len(mut_tar) > 0:
            ld_res.append(mut_bl)

    #record indepent networks with no mutal targets
    else:
        ind_net_lab=' + '.join(lds)
        ld_res.append(ind_net_lab)

#select edges with strongest p2p5 and remove duplicate and self edges
net_res_sort=net_res.sort_values('p2p5', ascending=False).drop('number_of_targets', axis=1)
net_res_filt=net_res_sort.loc[net_res_sort['Protein_name_A'] != net_res_sort['Protein_name_B']].drop_duplicates(['Protein_name_A', 'Protein_name_B'])

#recompute number of network targets for resolved networks and readjust to threshold
count_a=net_res_filt['S_ID_A'].value_counts()
count_df=pd.DataFrame({'S_ID_A': count_a.index, 'number_of_targets': count_a})
net_res_count=pd.merge(net_res_filt, count_df, on='S_ID_A', how='inner').sort_values('number_of_targets', ascending=False)
net_res_count=net_res_count.loc[net_res_count['number_of_targets'] >= 10]

#add resolved LD blocks to network summary
ld_res_df=pd.DataFrame({'LD_block': ld_un, 'LD_resolved': ld_res})
ld_pro_res=pd.merge(ld_pro, ld_res_df, on='LD_block', how='inner')

#export resolved network and resolved LD blocks
net_res_count.to_csv('results/findr/filtered/LD_resolution/AGES_findr_networks_adj_1FDR_net10_LD_resolved.tsv', sep='\t', index=None)
ld_pro_res.to_csv('results/findr/filtered/LD_resolution/AGES_findr_networks_adj_1FDR_net10_LD_resolved_summary.tsv', sep='\t', index=None)
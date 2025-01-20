import pandas as pd
import yaml

#read config
config_path='C:/INTRePID/scripts/config/network_analysis/ranked_networks.yml'

with open(config_path, 'r') as config_file:
    config=yaml.safe_load(config_file)

#import networks
net_fn=config['paths']['net_fn']
net=pd.read_csv(net_fn, delimiter='\t')

#import ranking and protein annotation
rank_fn=config['paths']['rank_fn']
rank=pd.read_csv(rank_fn, delimiter='\t')

def filt_net(n):

    #select top networks
    top=rank.loc[rank['Score'] >= n]

    #filter network
    net_top=pd.merge(net, top[['Protein_name_A']], on='Protein_name_A', how='inner')

    #count A and B-proteins
    a_count=net_top['Protein_name_A'].value_counts().to_dict()
    a_countdf=pd.DataFrame(a_count.items(), columns=['Protein_name_A', 'A-protein_count'])
    b_count=net_top['Protein_name_B'].value_counts().to_dict()
    b_countdf=pd.DataFrame(b_count.items(), columns=['Protein_name_B', 'B-protein_count'])

    #add counts to network
    net_a_count=pd.merge(net_top, a_countdf, on='Protein_name_A', how='left')
    net_ab_count=pd.merge(net_a_count, b_countdf, on='Protein_name_B', how='left')

    return net_ab_count

#filter and get A/B frequency
rank_n=config['filter']['rank']
net_rank=filt_net(rank_n)

#export ranked networks
out_pre=config['paths']['out_pre']
net_rank.to_csv(out_pre+'_ranked'+str(rank_n)+'.tsv', sep='\t', index=None)

#export ranked networks following B-protein filtering
if config['bfilt'] == True:
    btar=config['filter']['btar']
    net_bfilt=net_rank.loc[net_rank['B-protein_count'] >= btar]
    net_bfilt.to_csv(out_pre+'_ranked'+str(rank_n)+'_B'+str(btar)+'.tsv', sep='\t', index=None)
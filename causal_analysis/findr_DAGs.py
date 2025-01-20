import pandas as pd
import findr
import yaml

#read config
config_path='scripts/config/findr_DAGs.yml'

with open(config_path, 'r') as config_file:
        config = yaml.safe_load(config_file)

#import findr output and filtered networks
fin_fn=config['paths']['fin_fn']
fin=pd.read_csv(fin_fn, delimiter='\t')
net_fn=config['paths']['net_fn']
net=pd.read_csv(net_fn, delimiter='\t')

#select A-proteins from collapsed networks
ld_net=net[net['Protein_name_A'].str.contains(',')==True].drop_duplicates('Protein_name_A')

#get individual S_IDs from LD blocks
ld_sid=[]
for ld in ld_net['S_ID_A']:
	ld_split=(ld.split(','))
	for sid in ld_split:
		ld_sid.append(sid)

#get full list of somamers
full_sid=net.drop_duplicates('S_ID_A')['S_ID_A'].tolist()+ld_sid

#filter findr results for complete somamer list
fin_filt=fin.loc[(fin['S_ID_A'].isin(full_sid)) & (fin['S_ID_B'].isin(full_sid))]

#replace S_IDs with LD block labels
for index, row in ld_net.iterrows():
	sid_lab=row['S_ID_A']
	sid_spl=sid_lab.split(',')
	fin_filt=fin_filt.replace(to_replace=sid_spl, value=sid_lab)

#select edges with strongest p2p5, remove duplicate edges and set self edges to 0
fin_filt_sort=fin_filt.sort_values('p2p5', ascending=False).drop_duplicates(['S_ID_A', 'S_ID_B'])
fin_filt_sort.loc[fin_filt_sort['S_ID_A'] == fin_filt_sort['S_ID_B'], 'p2p5'] = 0

#convert from long to wide
fin_mat=fin_filt_sort.pivot(index='S_ID_A', columns='S_ID_B', values='p2p5')

#use findr to reconstruct acyclic graphs
l=findr.lib(path="/home/sean/temp/findr/libfindr.so",loglv=6,rs=0,nth=0)
ans=l.netr_one_greedy(fin_mat.values)

#convert from wide to long
ans_df=pd.DataFrame(ans['net'], index=fin_mat.index, columns=fin_mat.columns)
dag=ans_df.stack().reset_index()
dag.columns=['S_ID_A', 'S_ID_B', 'edge']

#remove self edges
dag_se=dag.loc[dag['S_ID_A'] != dag['S_ID_B']]

#select true edges
dag_true=dag_se.loc[dag_se['edge'] == True]

#filter for significant edges from original network
net_dag=pd.merge(net, dag_true, on=['S_ID_A', 'S_ID_B'], how='inner').drop(['edge', 'number_of_targets'], axis=1)

#add in B-protein edges to DAG network
aprot=list(set(net['Protein_name_A']))
net_bp=net.loc[~net['Protein_name_B'].isin(aprot)].drop('number_of_targets', axis=1)
dag_full=pd.concat([net_dag, net_bp])

#export results
out_fn=config['paths']['out_fn']
dag_full.to_csv(out_fn, sep='\t', index=None)
import pandas as pd
import numpy as np
import yaml

#read config
config_path="C:/INTRePID/scripts/config/network_analysis/networks_FDR_filter.yml"

with open(config_path, "r") as config_file:
        config=yaml.safe_load(config_file)

#import findr networks
res_fn=config['paths']['res_fn']
res=pd.read_csv(res_fn, delimiter='\t')
res_pp=res.loc[res['p2p5'] > 0.5]

#import protein reference
pro_ref=pd.read_csv('data/proteins/protein_annotations_filtered.txt', delimiter='\t')

#annotate with entrez ID
res_ref=pd.merge(res_pp, pro_ref[['S_ID', 'Protein_name']], left_on='S_ID_A', right_on='S_ID', how='left')
res_ref=pd.merge(res_ref, pro_ref[['S_ID', 'Protein_name']], left_on='S_ID_B', right_on='S_ID', how='left')

#reorientate dataframe
res_anno=res_ref[['S_ID_A', 'Target_A', 'Protein_name_x', 'S_ID_B', 'Target_B', 'Protein_name_y', 
        'p2', 'p3', 'p4', 'p5', 'p2p3', 'p2p5']]
res_anno.columns=['S_ID_A', 'Target_A', 'Protein_name_A', 'S_ID_B', 'Target_B', 'Protein_name_B', 
        'p2', 'p3', 'p4', 'p5', 'p2p3', 'p2p5']

#filter to remove non-human or QC somamers
res_filt=pd.merge(res_anno, pro_ref[['S_ID']], left_on='S_ID_A', right_on='S_ID', how='inner').drop('S_ID', axis=1)
res_filt=pd.merge(res_filt, pro_ref[['S_ID']], left_on='S_ID_B', right_on='S_ID', how='inner').drop('S_ID', axis=1)

#filter networks to selected FDR threshold
def filt_net(net, trd, test):
    steps=np.arange(0.5, 1, 0.01)
    fdr_l=[]
    for s in steps:
        fdr=1-net.loc[net[test] >= s][test].mean()
        fdr_l.append(fdr)
    fdr_filt=pd.DataFrame({test: steps, 'FDR': fdr_l})
    sel=fdr_filt.iloc[(fdr_filt['FDR']-trd).abs().argsort()[:1]].values[0][0]
    selr=np.round(sel, 2)
    net_fdr=net.loc[net[test] >= selr].sort_values(test, ascending=False)
    
    return net_fdr

fdr1=filt_net(res_filt, 0.01, 'p2p5')
fdr5=filt_net(res_filt, 0.05, 'p2p5')
fdr10=filt_net(res_filt, 0.1, 'p2p5')

#export FDR filtered networks
out_pre=config['paths']['out_pre']
fdr1.to_csv(out_pre+'1FDR.tsv', sep='\t', index=None)
fdr5.to_csv(out_pre+'5FDR.tsv', sep='\t', index=None)
fdr10.to_csv(out_pre+'10FDR.tsv', sep='\t', index=None)
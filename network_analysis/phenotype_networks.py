import pandas as pd
from multipy import fdr
import yaml

#read config
config_path="scripts/config/network_analysis/phenotype_networks.yml"

with open(config_path, "r") as config_file:
        config=yaml.safe_load(config_file)

#import network
net_fn=config['paths']['net_fn']
net=pd.read_csv(net_fn, delimiter='\t')

#import phenotype associations
ph_fn=config['paths']['ph_fn']
ph=pd.read_csv(ph_fn, delimiter='\t')

#adjust for multiple testing
qv=fdr.qvalue(ph['P'].values)
ph['q-value']=qv[1]

#filter for 5% FDR
qv=config['filter']['fdr']
ph_sig=ph.loc[ph['q-value'] <= qv]
ph_genes=ph_sig['EntrezGeneSymbol'].tolist()

#filter regulators as targets for phenotype associated genes
net_ph=net.loc[net['Protein_name_A'].isin(ph_genes) & net['Protein_name_B'].isin(ph_genes)]


#filter to remove duplicate somamers and self edges
un=config['filter']['unique']
if un == True:
    net_ph_filt=net_ph[net_ph.columns[:6]]
    net_ph_se=net_ph_filt.loc[net_ph_filt['Protein_name_A'] != net_ph_filt['Protein_name_B']]
    net_ph_out=net_ph_se.drop_duplicates(['Protein_name_A', 'Protein_name_B'])
else:
    net_reg_ph_out=net_ph[net_ph.columns[:6]]

#output results
out_fn=config['paths']['out_fn']
net_ph_out.to_csv(out_fn, sep='\t', index=None)
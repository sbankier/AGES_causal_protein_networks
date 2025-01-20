import pandas as pd
import yaml

#read config
config_path="scripts/config/network_analysis/regulators_as_targets.yml"

with open(config_path, "r") as config_file:
        config=yaml.safe_load(config_file)

#import networks
net_fn=config['paths']['net_fn']
fdr=pd.read_csv(net_fn, delimiter='\t')

#remove self-edges
fdr_filt=fdr[fdr['Protein_name_A'] != fdr['Protein_name_B']]

#filter network targets for other regulators
pro_a=fdr_filt[['Protein_name_A']].drop_duplicates()['Protein_name_A'].to_list()
regs=pd.DataFrame({'Protein_name_B': pro_a})
net_reg=pd.merge(fdr_filt, regs, on='Protein_name_B', how='inner')

#filter for unqiue proteins
net_reg_un=net_reg.drop_duplicates(['Protein_name_A', 'Protein_name_B'])

#output filtered network
out_fn=config['paths']['out_pre']+'_regs_as_tars.tsv'
net_reg_un.to_csv(out_fn, sep='\t', index=None)

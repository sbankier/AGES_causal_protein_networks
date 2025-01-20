import pandas as pd
import yaml

#read config
config_path='scripts/config/data_processing/protein_subsets.yml'

with open(config_path, 'r') as config_file:
    config=yaml.safe_load(config_file)

#import protein data
pro=pd.read_csv(config['paths']['input_fn'], delimiter='\t')

#sample parameters
n_sub=config['subsets']['n_sub']
size=config['subsets']['size']

#randomly sample individuals from protein data
rand_pro={}
for x in range(n_sub):
    sam=pro.sample(n=size, random_state=x)
    lab='sample_'+str(x+1)
    rand_pro[lab]=sam

#export subset dataframes
out_pre=config['paths']['output_pre']
for key, df in rand_pro.items():
    df.to_csv(out_pre+key+'.txt', sep='\t', index=False)
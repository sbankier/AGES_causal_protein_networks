import pandas as pd
import yaml

#read config
config_path='C:/INTRePID/scripts/config/network_analysis/network_target_filtering.yml'

with open(config_path, 'r') as config_file:
    config=yaml.safe_load(config_file)

#set target threshold
ntar=config['target_threshold']['ntar']

#import findr networks
net_pre=config['paths']['net_pre']
fdr1=pd.read_csv(net_pre+'1FDR.tsv', delimiter='\t')
fdr5=pd.read_csv(net_pre+'5FDR.tsv', delimiter='\t')
fdr10=pd.read_csv(net_pre+'10FDR.tsv', delimiter='\t')

#organise networks by size
def count_targets(fdr):
    count_a=fdr['S_ID_A'].value_counts()
    count_df=pd.DataFrame({'S_ID_A': count_a.index, 'number_of_targets': count_a})
    fdr_count=pd.merge(fdr, count_df, on='S_ID_A', how='inner').sort_values('number_of_targets', ascending=False)
    return fdr_count

fdr1_count=count_targets(fdr1)
fdr5_count=count_targets(fdr5)
fdr10_count=count_targets(fdr10)

#filter networks based on number of targets
fdr1_filt=fdr1_count.loc[fdr1_count['number_of_targets'] >= ntar]
fdr5_filt=fdr5_count.loc[fdr5_count['number_of_targets'] >= ntar]
fdr10_filt=fdr10_count.loc[fdr10_count['number_of_targets'] >= ntar]

#output filtered networks
out_pre=config['paths']['out_pre']
fdr1_filt.to_csv(out_pre+'1FDR_net'+str(ntar)+'.tsv', sep='\t', index=None)
fdr5_filt.to_csv(out_pre+'5FDR_net'+str(ntar)+'.tsv', sep='\t', index=None)
fdr10_filt.to_csv(out_pre+'10FDR_net'+str(ntar)+'.tsv', sep='\t', index=None)
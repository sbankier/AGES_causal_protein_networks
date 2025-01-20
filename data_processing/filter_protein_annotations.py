import pandas as pd
import numpy as np

#import annotation file and filter for human proteins/ mouse fusion and remove QC probes
ano=pd.read_csv('data/proteins/protein_annotations_base.txt', delimiter='\t')
ano_filt=ano.loc[ano['PROBE_CLASS'] == 'HUMAN_SELEX']

#obtain same SEQ_ID annotations as in 7K soma data
new_sid=[]
for x in ano_filt['SeqId']:
    splt=x.split('_')
    sid='seq.'+splt[0]+'.'+splt[1]
    new_sid.append(sid)

#reorientate columns
ano_filt['S_ID']=new_sid
cols=ano_filt.columns.tolist()
cols=cols[-1:]+cols[:-1]
ano_filt=ano_filt[cols]

#replace "None" with NaN
ano_filt['EntrezGeneSymbol'] = ano_filt['EntrezGeneSymbol'].replace('None', np.nan)

#include additional annotation column where if "EntrezGeneSymbol" is unavailable use "SYMBOL"
new_col=ano_filt['EntrezGeneSymbol'].fillna(ano_filt['SYMBOL'])
ano_filt['Protein_name']=new_col

#export filtered annotation file
ano_filt.to_csv('data/proteins/protein_annotations_filtered.txt', sep='\t', index=None)

#import protein data
pro_base=pd.read_csv('data/proteins/xpsc_a1_full_nonhuman.txt', delimiter='\t')
pro_fu=pd.read_csv('data/proteins/xpsc_a2_full_nonhuman.txt', delimiter='\t')

#filter protein measurments for filtered annotations
def filter_pro(pro):
    filt_sid=['S_ID', 'PROJ']+ano_filt['S_ID'].tolist()
    pro_filt=pro[filt_sid]
    return pro_filt

pro_base_filt=filter_pro(pro_base)
pro_fu_filt=filter_pro(pro_fu)

#export filtered protein data
pro_base_filt.to_csv('data/proteins/xpsc_a1_full.txt', sep='\t', index=None)
pro_fu_filt.to_csv('data/proteins/xpsc_a2_full.txt', sep='\t', index=None)
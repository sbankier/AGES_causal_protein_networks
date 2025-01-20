import os
import numpy as np
import pandas as pd
import pickle
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score

#import full network matrix of all possible edges
full=pd.read_csv('results/findr/AGES_cispQTLs_5FDR_findr_output_adj-age-sex.tsv', delimiter='\t')
fullp=full.pivot(index='S_ID_A', columns='S_ID_B', values='p2p5')
full_mat=pd.DataFrame(0, index=fullp.index, columns=fullp.columns)

#obtain a flattened matrix of true edges from a given network
def get_edges(df):
    pred_mat=full_mat.copy()
    for index, row in df.iterrows():
        sida=row['S_ID_A']
        sidb=row['S_ID_B']
        pred_mat.at[sida, sidb]=1
    pred_out=pred_mat.values.flatten()
    
    return pred_out

#import ground truth networks
print('\n#######################\n# IMPORT GROUND TRUTH #\n#######################\n')
gt_edges_fdr={}
for f in ['1FDR', '5FDR', '10FDR']:
    net=pd.read_csv('results/findr/FDR/AGES_findr_networks_adj_'+f+'.tsv', delimiter='\t')
    edges=get_edges(net)
    gt_edges_fdr[f]=edges

#get all sub-sampled network paths
def list_files(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
sub_path='results/findr/random_sampling/'
sub_dir=list_files(sub_path)

#get findr scores from sub-sampled networks
print('\n############################\n# GET SUB-SAMPLED NETWORKS #\n############################\n')
mat_set={}
for f in sub_dir:
    #get findr scores for sub-sampled networks as flattened matrix
    net=pd.read_csv(sub_path+f, delimiter='\t')
    netp=net.pivot(index='S_ID_A', columns='S_ID_B', values='p2p5')
    net_mat=netp.values.flatten()
    net_mat[np.isnan(net_mat)]=0
    mat_set.update({f: net_mat})
    print(f)

#obtain sets of test and ground truth scores
sub_set={}
for x in mat_set:
    fn_sp=x.split('_')
    n=fn_sp[7]
    sam=fn_sp[9].split('.')[0]
    sub_set.update({x+'1FDR': ('1FDR', n, sam,  mat_set[x], gt_edges_fdr['1FDR'])})
    sub_set.update({x+'5FDR': ('5FDR', n, sam, mat_set[x], gt_edges_fdr['5FDR'])})
    sub_set.update({x+'10FDR': ('10FDR', n, sam, mat_set[x], gt_edges_fdr['10FDR'])})

#calculate ROC AUC for sub-sampled test sets
print('\n#################\n# CALCULATE AUC #\n#################\n')
n_l=[]
fdr_l=[]
sam_l=[]
auc_l=[]
for fn, (fdr, n, sam, y_test, y_true) in sub_set.items():
    n_l.append(n)
    fdr_l.append(fdr)
    sam_l.append(sam)
    auc=roc_auc_score(y_true=y_true, y_score=y_test)
    auc_l.append(auc)
    print(fn)

#assemble ROC results as dataframe
auc_df=pd.DataFrame({'n': n_l, 'sample': sam_l, 'FDR': fdr_l, 'AUC': auc_l})

#calculate AUC mean and std deviation
auc_mean_l=[]
auc_std_l=[]
for index, row in auc_df.iterrows():
    auc_filt=auc_df.loc[(auc_df['n'] == row['n']) & (auc_df['FDR'] == row['FDR'])]
    auc_mean=auc_filt['AUC'].mean()
    auc_std=auc_filt['AUC'].std()
    auc_mean_l.append(auc_mean)
    auc_std_l.append(auc_std)
auc_df['AUC_mean']=auc_mean_l
auc_df['AUC_std']=auc_std_l

#export AUC results
auc_df.to_csv('results/findr/random_sampling/results/AGES_findr_random_FDR_AUC.tsv', sep='\t', index=None)

#calculate precision recall for test sets
print('\n##############################\n# CALCULATE PRECISION-RECALL #\n##############################\n')
pr_set={}
for fn, (fdr, n, sam, y_test, y_true) in sub_set.items():
    p, r, t = precision_recall_curve(y_true=y_true, y_score=y_test, drop_intermediate=True)
    a = average_precision_score(y_true=y_true, y_score=y_test)
    pr_set.update({fn: (fdr, n, sam, p, r, t, a)})
    print(fn)

#export precision recall results
with open('results/findr/random_sampling/results/AGES_findr_random_FDR_PR.pkl', 'wb') as f:
    pickle.dump(pr_set, f)
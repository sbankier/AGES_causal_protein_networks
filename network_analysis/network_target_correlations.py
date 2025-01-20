import pandas as pd
import numpy as np
import math
from scipy import stats
from multipy.fdr import qvalue
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

#import FDR and target filtered protein networks and get list of network A-proteins
fdr1=pd.read_csv('results/findr/filtered/AGES_findr_networks_adj_1FDR_net10.tsv', delimiter='\t')
net_list=fdr1[['S_ID_A']].drop_duplicates()['S_ID_A'].values

#function to transpose dataframe
def flip_df(df):
    df=df.T.reset_index()
    df.columns = df.iloc[0]
    df=df.iloc[1:]
    return df

#import protein expression data
pro=pd.read_csv('data/proteins/xpsc_a1_full_adj-sex-age.txt', delimiter='\t').dropna()
pro_mat=flip_df(pro.drop('PROJ', axis=1))

#for all protein networks perform Kruskal Wallis test for network target vs random pairwise correlations
kw_l=[]
for x in net_list:
    #get absolute correlation values as matrix
    net=fdr1.loc[fdr1['S_ID_A'] == x]
    net_pro=pd.merge(net[['S_ID_B']], pro_mat, left_on='S_ID_B', right_on='S_ID').drop('S_ID', axis=1)
    net_mat=flip_df(net_pro)
    corr=net_mat.drop('S_ID_B', axis=1).astype(float).corr().abs().reset_index().fillna(0)
    corr.rename(columns={0:'S_ID_B'}, inplace=True)

    #get random correlation values as matrix
    rn=len(net)
    net_pro_rand=pro_mat.sample(n=rn, random_state=100)
    net_mat_rand=flip_df(net_pro_rand)
    corr_rand=net_mat_rand.drop('S_ID', axis=1).astype(float).corr().abs().reset_index().fillna(0)
    corr_rand.rename(columns={0:'S_ID_B'}, inplace=True)

    #remove diagonal from correlation matrices
    def rmv_diag(mat):
        matnan=mat.replace(1.0,np.NaN).drop('S_ID_B', axis=1).values
        mat_nodiag=matnan[np.logical_not(np.isnan(matnan))]
        return mat_nodiag

    corr_nodiag=rmv_diag(corr)
    corr_rand_nodiag=rmv_diag(corr_rand)

    #conduct Kruskal Wallis test between target and random correlations
    kw=stats.kruskal(corr_nodiag, corr_rand_nodiag)
    kw_l.append(kw)
    print(x)

#format Kruskal Wallis results as dataframe
kwdf=pd.DataFrame(kw_l)
kwdf['S_ID_A']=net_list

#correct for multiple testing
qv=qvalue(kwdf['pvalue'].values)
kwdf['q-value']=qv[1]
bf=multipletests(kwdf['pvalue'].values, method='bonferroni')
kwdf['bf']=bf[1]

#re-order columns
kwdf_out=kwdf[['S_ID_A', 'statistic', 'pvalue', 'q-value', 'bf']]
kwdf_out.columns=['S_ID_A', 'KW_statistic', 'p-value', 'q-value', 'bonferroni']

#export KW test results
kwdf_out.to_csv('results/findr/network_correlations/AGES_findr_networks_1FDR_netvrand_KW.tsv', sep='\t', index=None)
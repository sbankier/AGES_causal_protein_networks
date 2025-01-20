import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#import FDR and target filtered protein networks and get list of network A-proteins
fdr1=pd.read_csv('C:/INTRePID/results/findr/filtered/AGES_findr_networks_adj_1FDR_net10.tsv', delimiter='\t')
net_list=fdr1[['S_ID_A']].drop_duplicates()['S_ID_A'].values

#import protein expression data
pro=pd.read_csv('data/proteins/xpsc_a1_full_adj-sex-age.txt', delimiter='\t').dropna()

#function for obtaining eignevectors from exisiting protein networks with randomised targets
def get_ev(v):

    #create dictionary to store results
    pcdf={"pc1df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "pc2df": pd.DataFrame({'AGES_ID': pro['S_ID']}), 
    "pc3df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "pc4df": pd.DataFrame({'AGES_ID': pro['S_ID']}), 
    "pc5df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "var_1": [], "var_2": [], "var_3": [], "var_4": [], "var_5": []}

    #obtain first 5 principal components and corresponding variance explained for each sub-network
    for x in net_list:
        #obtain protein measurements for network proteins
        net_len=len(fdr1.loc[fdr1['S_ID_A'] == x])
        pro_b=pd.DataFrame({"S_ID_B": pro.columns[2:]})
        net_rand=pro_b.sample(net_len)['S_ID_B'].tolist()
        net_rand.append(x)
        net_pro=pro[net_rand]
        net_pro

        #apply standard scaler to dataset features
        ct_x = net_pro.values
        ct_xt = StandardScaler().fit_transform(ct_x)

        #Obtain first 5 prinicipal components
        ct_pca = PCA(n_components=5)
        ct_principalComponents = ct_pca.fit_transform(ct_x)
        pcdf['pc1df'][x]=ct_principalComponents[:,0]
        pcdf['pc2df'][x]=ct_principalComponents[:,1]
        pcdf['pc3df'][x]=ct_principalComponents[:,2]
        pcdf['pc4df'][x]=ct_principalComponents[:,3]
        pcdf['pc5df'][x]=ct_principalComponents[:,4]

        #get variance explained by each PC
        var_exp=ct_pca.explained_variance_ratio_
        pcdf['var_1'].append(round(var_exp[0]*100, 2))
        pcdf['var_2'].append(round(var_exp[1]*100, 2))
        pcdf['var_3'].append(round(var_exp[2]*100, 2))
        pcdf['var_4'].append(round(var_exp[3]*100, 2))
        pcdf['var_5'].append(round(var_exp[4]*100, 2))


    #format variance explained as a dataframe
    vardf=pd.DataFrame({'S_ID_A': net_list, 'PC1_variance_explained': pcdf['var_1'], 'PC2_variance_explained': pcdf['var_2'], 'PC3_variance_explained': pcdf['var_3'],
    'PC4_variance_explained': pcdf['var_4'], 'PC5_variance_explained': pcdf['var_5']})

    #export output
    pcdf['pc1df'].to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_eigenvectors_PC1_'+v+'.tsv', sep='\t', index=None)
    pcdf['pc2df'].to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_eigenvectors_PC2_'+v+'.tsv', sep='\t', index=None)
    pcdf['pc3df'].to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_eigenvectors_PC3_'+v+'.tsv', sep='\t', index=None)
    pcdf['pc4df'].to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_eigenvectors_PC4_'+v+'.tsv', sep='\t', index=None)
    pcdf['pc5df'].to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_eigenvectors_PC5_'+v+'.tsv', sep='\t', index=None)
    vardf.to_csv('C:/INTRePID/results/findr/random_network_eigenvectors/AGES_findr_random_networks_adj_1FDR_net10_PC1-5_var-exp_'+v+'.tsv', sep='\t', index=None)

#run for multiple iterations
iter_v=['v1', 'v2', 'v3', 'v4', 'v5']

for v in iter_v:
    get_ev(v)
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
import yaml

#read config
config_path="C:/INTRePID/scripts/config/network_analysis/networks_eigenvectors.yml"

with open(config_path, "r") as config_file:
    config=yaml.safe_load(config_file)

#import FDR and target filtered protein networks and get list of network A-proteins
res_fn=config['paths']['res_fn']
fdr1=pd.read_csv(res_fn, delimiter='\t')
net_list=fdr1[['S_ID_A']].drop_duplicates()['S_ID_A'].values

#import protein expression data
pro_fn=config['paths']['pro_fn']
pro=pd.read_csv(pro_fn, delimiter='\t').dropna()

#create dictionary to store results
pcdf={"pc1df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "pc2df": pd.DataFrame({'AGES_ID': pro['S_ID']}), 
"pc3df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "pc4df": pd.DataFrame({'AGES_ID': pro['S_ID']}), 
"pc5df": pd.DataFrame({'AGES_ID': pro['S_ID']}), "var_1": [], "var_2": [], "var_3": [], "var_4": [], "var_5": []}

#obtain first 5 principal components and corresponding variance explained for each sub-network
for x in net_list:

    #obtain protein measurements for network proteins
    net=fdr1.loc[fdr1['S_ID_A'] == x]['S_ID_B'].to_list()
    net.append(x)

    #separate combined S_IDs
    sids=[]
    for sid in net:
        if ',' in sid:
            sids.extend(sid.split(','))
        else:
            sids.append(sid)

    net_pro=pro[sids]

    #apply standard scaler to dataset features
    ct_x=net_pro.values
    ct_xt=StandardScaler().fit_transform(ct_x)

    #Obtain first 5 prinicipal components
    ct_pca=PCA(n_components=5)
    ct_principalComponents=ct_pca.fit_transform(ct_x)
    pc1=ct_principalComponents[:,0]
    pc2=ct_principalComponents[:,1]
    pc3=ct_principalComponents[:,2]
    pc4=ct_principalComponents[:,3]
    pc5=ct_principalComponents[:,4]

    #format sign of the PCs depending on correlation with A-protein
    def align_sign(pc, aprot):
        
        #if unresolved then take mean expression for correlation
        if ',' in aprot:
            ap=aprot.split(',')
            ap_exp=pro[ap].mean(axis=1).values
        else:
            ap_exp=pro[aprot].values

        #alternate sign if correlation is negative    
        cor=stats.spearmanr(pc, ap_exp)
        if cor[0] > 0:
            n_pc=pc
        elif cor[0] < 0:
            n_pc=-pc

        return n_pc

    #export PCs
    pcdf['pc1df'][x]=align_sign(pc1, x)
    pcdf['pc2df'][x]=align_sign(pc2, x)
    pcdf['pc3df'][x]=align_sign(pc3, x)
    pcdf['pc4df'][x]=align_sign(pc4, x)
    pcdf['pc5df'][x]=align_sign(pc5, x)

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

#annotate with protein names
ano=fdr1[['S_ID_A', 'Protein_name_A']].drop_duplicates()
vardf_out=pd.merge(ano, vardf, on='S_ID_A', how='inner')

#export output
out_pre=config['paths']['out_pre']
pcdf['pc1df'].to_csv(out_pre+'_eigenvectors_PC1.tsv', sep='\t', index=None)
pcdf['pc2df'].to_csv(out_pre+'_eigenvectors_PC2.tsv', sep='\t', index=None)
pcdf['pc3df'].to_csv(out_pre+'_eigenvectors_PC3.tsv', sep='\t', index=None)
pcdf['pc4df'].to_csv(out_pre+'_eigenvectors_PC4.tsv', sep='\t', index=None)
pcdf['pc5df'].to_csv(out_pre+'_eigenvectors_PC5.tsv', sep='\t', index=None)
vardf_out.to_csv(out_pre+'_PC1-5_var-exp.tsv', sep='\t', index=None)
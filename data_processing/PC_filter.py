import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#import protein data
a1_matched=pd.read_csv('data/proteins/xpsc_a1_a2-matched_adj-sex-age.txt', delimiter='\t')
a2_matched=pd.read_csv('data/proteins/xpsc_a2_a1-matched_adj-sex-age.txt', delimiter='\t')

#conduct PCA
def get_pca(ct):
    #apply standard scaler to dataset features
    ct_x = ct.drop(['S_ID', 'PROJ'], axis=1).dropna().values
    ct_xt = StandardScaler().fit_transform(ct_x)

    #reduce data across 2 principal components
    ct_pca = PCA(n_components=2)
    ct_principalComponents = ct_pca.fit_transform(ct_x)
    ct_principalDf = pd.DataFrame(data = ct_principalComponents
                 , columns = ['principal component 1', 'principal component 2'])

    return ct_principalDf

a2_pc=get_pca(a2_matched)

#Find and remove outliers
pcdf=pd.concat([a2_matched[['S_ID']], a2_pc], axis=1)
pc_sort=pcdf.sort_values('principal component 2', ascending=False)
pc_filt=pc_sort.iloc[1:]

#remove outliers from original protein data
a1_filt=pd.merge(a1_matched, pc_filt[['S_ID']], on='S_ID', how='inner')
a2_filt=pd.merge(a2_matched, pc_filt[['S_ID']], on='S_ID', how='inner')

#export
a1_filt.to_csv('data/proteins/xpsc_a1_a2-matched_filt_adj-sex-age.txt', sep='\t', index=None)
a2_filt.to_csv('data/proteins/xpsc_a2_a1-matched_filt_adj-sex-age.txt', sep='\t', index=None)
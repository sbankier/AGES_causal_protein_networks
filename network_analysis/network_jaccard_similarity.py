import pandas as pd
from itertools import combinations
from tqdm import tqdm

#import network and protein_annotations
net=pd.read_csv('results/findr/filtered/LD_resolution/AGES_findr_networks_adj_1FDR_net10_LD_resolved.tsv', delimiter='\t')
net_un=net.drop_duplicates(['Protein_name_A', 'Protein_name_B'])

#obtain jaccord coefficent of similarity between sub-network targets
def get_jaccard(a1, a2):

    #select sub-networks
    sub1=net_un.loc[net_un['Protein_name_A'] == a1]
    sub2=net_un.loc[net_un['Protein_name_A'] == a2]

    #get representation of network within union
    sub1b=set(sub1['Protein_name_B'].tolist())
    sub2b=set(sub2['Protein_name_B'].tolist())

    #get union and intersection of targets
    union=sub1b.union(sub2b)
    intrsec=sub1b.intersection(sub2b)

    #calculate jaccard score
    js=len(intrsec)/len(union)

    return js

#generate all pairs of A-proteins
aprot=net_un['Protein_name_A'].unique()
aprot_pairs = list(combinations(aprot, 2))

#obtain jaccard score for all sub-network pairs
js_l=[]
for pair in tqdm(aprot_pairs):
    a1, a2 = pair
    js=get_jaccard(a1, a2)
    js_l.append(js)

#assmble output as dataframe
jsdf=pd.DataFrame(aprot_pairs, columns=['A-protein_1','A-protein_2'])
jsdf['jaccard_score']=js_l

#export results
jsdf.to_csv('results/findr/filtered/subnetwork_similarity/AGES_findr_networks_adj_1FDR_net10_LD_resolved_jaccard.tsv', sep='\t', index=None)
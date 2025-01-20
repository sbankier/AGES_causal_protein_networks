import pandas as pd
import yaml

#read config
config_path='scripts/config/network_analysis/subnetwork_summary.yml'

with open(config_path, "r") as config_file:
        config=yaml.safe_load(config_file)

#improt parameters
fdr=config['parameters']['FDR']
pro=config['parameters']['pro']

#report title
print('\n\n### '+pro+' '+str(fdr)+'%'+' FDR sub-network summary'+' ###\n')

#import findr networks
nets=pd.read_csv('C:/INTRePID/results/findr/FDR/AGES_findr_networks_adj_'+str(fdr)+'FDR.tsv', delimiter='\t')
net=nets.loc[nets['Protein_name_A'] == pro]

#organise networks by size
count_a=net['S_ID_A'].value_counts()
count_df=pd.DataFrame({'S_ID_A': count_a.index, 'number_of_targets': count_a})
net_count=pd.merge(count_df, net, on='S_ID_A', how='inner')

#report number of somamers
sid=net_count.sort_values('number_of_targets', ascending=False).drop_duplicates('S_ID_A')
sid_top=sid['S_ID_A'].values[0]
sid_n=len(sid)

if sid_n == 1:
    print(pro+' has 1 sommaer at '+str(fdr)+'%'+' FDR: '+sid_top)
else:
    print(pro+' has '+str(sid_n)+' sommaers at '+str(fdr)+'%'+' FDR. Selecting largest network driven by '+sid_top)
    net=net.loc[net['S_ID_A'] == sid_top]

#report number of targets
edge_n=len(net.drop_duplicates('Protein_name_B'))
print('number of unique targets: '+str(edge_n)+'\n')

#import network associations and reference
aprot=pd.read_csv('results/pheno_assoc/A-proteins/Aprot_assoc_5FDR_intersection.tsv', delimiter='\t')
net_ap=aprot.loc[aprot['SEQ_ID'] == sid_top]
eign=pd.read_csv('results/pheno_assoc/PC1_var8_assoc_5FDR_intersection.tsv', delimiter='\t')
net_ei=eign.loc[eign['xpcol_seqid'] == sid_top]

#filter for significant associations
ap_ph=net_ap.columns[net_ap.iloc[0] == 1]
ei_ph=net_ei.columns[net_ei.iloc[0] == 1]

#report A-protein associations
if len(ap_ph) == 0:
    print(pro+' is not associated with any AGES phenotype\n')
else:
    print(pro+' is associated with '+str(len(ap_ph))+' AGES phenotypes at a 5'+'%'+ ' FDR threshold:')
    for x in ap_ph:
        ph=x.split('_')[0]
        print('\t'+ph)

#report eigenprotein assocations
if len(ei_ph) == 0:
    print('\n'+pro+' network eigenprotein is not associated with any AGES phenotype\n')
else:
    print('\n'+pro+' network eigenprotein is associated with '+str(len(ei_ph))+' AGES phenotypes at a 5'+'%'+ 'FDR threshold:')
    for x in ei_ph:
        ph=x.split('_')[0]
        print('\t'+ph)

#import functional enrichment results
go=pd.read_csv('C:/INTRePID/results/functional_enrichment/networks_adj_'+str(fdr)+'FDR_net10_GO_terms.tsv', delimiter='\t')
go_net=go.loc[go['S_ID_A'] == sid_top].sort_values('p_value')
go_n=len(go_net)

#report enrichment results
if go_n == 0:
    print('\n'+pro+' sub-network is not enriched for any terms')
elif go_n > 0 and go_n < 11:
    print('\n'+pro+' is enriched for '+str(go_n)+' terms:')
    for x in go_net['term_name']:
        print('\t'+x)
else:
    print('\n'+pro+' sub-network is enriched for '+str(go_n)+ ' terms. The top 10 are:')
    for x in go_net.head(10)['term_name']:
        print('\t'+x)
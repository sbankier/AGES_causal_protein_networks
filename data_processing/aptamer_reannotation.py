import pandas as pd
import yaml

#read config
config_path='C:/INTRePID/scripts/config/data_processing/aptamer_reannotation.yml'

with open(config_path, 'r') as config_file:
    config=yaml.safe_load(config_file)

#import aptamer file
apt_fn=config['paths']['apt_fn']
apt=pd.read_csv(apt_fn, delimiter='\t').reset_index()

#select combined and single aptamers
apt_un=apt.loc[~apt['Aptamer'].str.contains(',')]
apt_com=apt.loc[apt['Aptamer'].str.contains(',')]

#re-annotate single aptamers
un_l=[]
for x in apt_un['Aptamer']:
    sp=x.split('-')
    n_lab='.'.join(['seq',sp[0],sp[1]])
    un_l.append(n_lab)

apt_un['S_ID']=un_l

#re annotate combined aptamers
def com_lab(com_str):
    aptamers=com_str.split(',')
    n_lab=[f"seq.{aptamer.replace('-', '.')}" for aptamer in aptamers]
    return ','.join(n_lab)

apt_com['S_ID']=apt_com['Aptamer'].apply(com_lab)

#combine results as single dataframe
n_apt=pd.concat([apt_un, apt_com]).sort_values('index').drop('index', axis=1)

#set S_ID as first column
sid=n_apt.pop('S_ID')
n_apt.insert(0, 'S_ID', sid)

#export results
out_pre=config['paths']['out_fn']
n_apt.to_csv(out_pre, sep='\t', index=None)
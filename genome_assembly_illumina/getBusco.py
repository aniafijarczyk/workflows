import glob
import pandas as pd
import re

buscos = glob.glob('run_busco_*/short_summary_busco_*.txt')
species = [ele.split('/')[0].replace('run_busco_','') for ele in buscos]
len(species)

singles = glob.glob('run_busco_*/single_copy_busco_sequences')
len(singles)
singles_sp = [ele.split('/')[0].replace('run_busco_','') for ele in singles]
singles_sp

R = []
for spec in species:
    fname = 'run_busco_'+spec+'/short_summary_busco_'+spec+'.txt'
    fh = open(fname,'r').readlines()
    k = [ele.split() for ele in fh if ele[0] not in '#'][1][0]
    kr = re.split(',|\[|\]',k)
    D = {a.split(':')[0]:[float(a.split(':')[1].replace('%',''))] for a in kr if a}
    D['species'] = [spec]
    df = pd.DataFrame.from_dict(D)
    R.append(df)
dR = pd.concat(R).reset_index().drop(columns=['index'])
dR.head()
dR.to_csv('getBusco.tab',sep='\t',header=True,index=False)

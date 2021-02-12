import pandas as pd
import glob
from Bio import SeqIO
import sys

#fname = 'blastx_TB5_uniprot_merged.tab'
fname = sys.argv[1]
#fai = 'spades_TB5_clean_softmasked.fasta.fai'
fai = sys.argv[2]


F = []
for fai in [fai]:
    #fai=fais[0]
    #name = fai.split("/")[-1].replace(".fna.fai","")
    name = fai.split(".")[0]
    dlen = pd.read_csv(fai, sep="\t", header=None, names=['chrom','len','x1','x2','x3'])

    df = pd.read_csv(fname,
                     sep="\t",header=None, names=['chrom','hit_name','id','length','mismatch','gapopen','qstart',
                                                  'qstop','sstart','sstop','evalue','bitscore','salltitles','taxid',
                                                  'start','stop','group','taxon','new_group'])
    df['hlen'] = df['stop'] - df['start']
    dg = df.groupby(['chrom','taxon']).agg({'hlen':'sum'}).reset_index()
    df1 = dg.loc[dg['taxon'] == 'fungi',['chrom','hlen']]
    df2 = dg.loc[dg['taxon'] == 'bacvir',['chrom','hlen']]
    dm1 = pd.merge(dlen,df1,on=['chrom'],how='left').fillna(0).rename(columns={'hlen':'fun_len'})
    dm2 = pd.merge(dm1,df2,on=['chrom'],how='left').fillna(0).rename(columns={'hlen':'bacvir_len'})
    dm2['fun_frac'] = dm2['fun_len']/dm2['len']
    dm2['bacvir_frac'] = dm2['bacvir_len']/dm2['len']
    dm2['species'] = name
    dm3 = dm2[['chrom','fun_frac','bacvir_frac','len','species']]
    F.append(dm3)
    
dF = pd.concat(F)
dF['cont_frac'] = dF['fun_frac'] - dF['bacvir_frac']
dF['cont'] = dF.apply(lambda x: 'contaminated' if x['fun_frac'] < x['bacvir_frac'] else 'ok', axis=1)
dsF = dF.sort_values(by=['cont_frac'],ascending=True)
dsF.to_csv('calculateOverlap.tab',sep='\t',header=True,index=False)


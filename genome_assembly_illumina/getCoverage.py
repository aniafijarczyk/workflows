import pandas as pd
import glob

deps = glob.glob("../depths/*.out")

N = []
for fname in deps:
    name = fname.split("/")[-1].replace("_depth.out","")
    print(name)
    cov = pd.read_csv(fname,sep='\t',header=None,names=['chrom','pos','cov'])
    medianCov = cov['cov'].median()
    #cov['ncov'] = cov['cov']/medianCov
    covsum = cov.groupby(['chrom']).agg({'cov':['count','mean','std']}).reset_index()
    covsum.columns = covsum.columns.map('_'.join)
    covsum['ncov'] = covsum['cov_mean']/medianCov
    covsum['species'] = name
    covsum['chrom'] = covsum['chrom_'].apply(lambda x: x.split("|")[0])
    novsum = covsum[['species','chrom','cov_mean','ncov']]
    N.append(novsum)
dN = pd.concat(N)

dN.to_csv('getCoverage.tab', header=True, index=False, sep='\t')

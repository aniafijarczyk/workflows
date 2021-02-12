import pandas as pd
import glob
import sys
from collections import defaultdict
from operator import itemgetter
import time
'''
The script reads merged formated output of blast with best hits. It finds groups of overlapping hits and for each group selects best hit (best bitscore) and outputs only best hits. Columns 'start' and 'stop' which are the beginning and end of best hit in the group, 'gropup' - number of the group, and 'taxon' which is the name of the taxon set in the file name.
'''

filename = sys.argv[1]

#tabs = glob.glob("blastx_*_merged_bestoverlap2.tab")

start_time = time.time()

for tab in [filename]:
    #tab = tabs[0]
    print(tab)
    #name = tab.split("/")[-1].replace("blastx_","").replace(".tab","")
    name = tab.split("/")[-1].split(".")[0]
    df = pd.read_csv(tab, sep="\t", header=None, names=['chrom','hit_name','id','length','mismatch',
                                                        'gapopen','qstart','qstop','sstart','sstop',
                                                        'evalue','bitscore','salltitles','taxid','start',
                                                        'stop','group','taxon'])
    dftab = df.values.tolist()
    d = defaultdict(list)
    for ele in dftab:
        d[ele[0]].append(ele)

    R = []
    for ref in d.keys():
        BigG = []
        # sort by start position
        dss = sorted(d[ref],key=itemgetter(14))

        #print('getting overlaps')
        D = set([])
        G = []
        i = 0
        for ele in dss:
            pos1 = ele[14]
            pos2 = ele[15]
            interv = set(range(pos1,pos2+1))
            if D.intersection(interv):
                G.append(ele+[i])
                D = D.union(interv)
            else:
                i+=1
                G.append(ele+[i])
                D = D.union(interv)
        BigG+=G

        B = []
        dgr = defaultdict(list)
        for hit in BigG:
            dgr[hit[-1]].append(hit)

        for gr in dgr.keys():
            # sort each group from highest to lowest bitscore & take the best
            dgsort = sorted(dgr[gr],key=itemgetter(11),reverse=True)
            dbest = dgsort[0]
            B.append(dbest)
        R+=B

    wh = open(str(name)+'_bestoverlap.tab','w')
    for ele in R:
        f = '\t'.join([str(j) for j in ele])+'\n'
        wh.write(f)
    wh.flush()
    wh.close()
fin = time.time() - start_time

    

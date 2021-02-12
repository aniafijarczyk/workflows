import pandas as pd
import glob
import sys
from collections import defaultdict
from operator import itemgetter
import time

'''
The script reads blastx output for blast of long contigs against protein database. It finds groups of overlapping hits and for each group selects best hit (best bitscore) and outputs only best hits. Additional columns include 'start' and 'stop' which are the beginning and end of cumulative set of hits from the group, 'gropup' - number of the group, and 'taxon' which is the name of the taxon set in the file name.
'''

code = sys.argv[1] # a name tag attacked to an output in last column
sampfile = sys.argv[2]

#tabs = glob.glob("blastx_*_uniprot_"+code+".tab")

for tab in [sampfile]:
    #tab = tabs[0]
    print(tab)
    #name = tab.split("/")[-1].replace("blastx_","").replace(".tab","")
    name = tab.split('/')[-1].split('.')[0]
    df = pd.read_csv(tab, sep="\t", header=None, names=['chrom','hit_name','id','length','mismatch','gapopen','qstart',
                                                           'qstop','sstart','sstop','evalue','bitscore','salltitles','taxid'])

    dftab = df.values.tolist()
    d = defaultdict(list)
    for ele in dftab:
        d[ele[0]].append(ele)

    R = []
    for ref in d.keys():
        BigG = []
        ds = []
        for hsp in d[ref]:
            start = min(hsp[6],hsp[7])
            stop = max(hsp[6],hsp[7])
            ds.append([start,stop] + hsp)
        dss = sorted(ds,key=itemgetter(0))

        #print('getting overlaps')
        D = set([])
        G = []
        i = 0
        for ele in dss:
            pos1 = ele[0]
            pos2 = ele[1]
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
            dgsort = sorted(dgr[gr],key=itemgetter(13),reverse=True)
            dbest = dgsort[0]
            # get min start and max stop for the group and add to best hit
            dgstart = min([ele[0] for ele in dgr[gr]])
            dgstop = max([ele[1] for ele in dgr[gr]])
            dbest = dbest[2:-1] + [dgstart,dgstop,dbest[-1],code]
            B.append(dbest)
        R+=B

    wh = open(str(name)+'_bestoverlap.tab','w')
    for ele in R:
        f = '\t'.join([str(j) for j in ele])+'\n'
        wh.write(f)
    wh.flush()
    wh.close()

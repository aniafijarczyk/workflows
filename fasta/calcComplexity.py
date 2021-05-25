from Bio import SeqIO
import re
from itertools import permutations
import sys

# The script estimates complexity for a genetic sequence of a given range around selected mutations
# The script reads: 
# 1) genome in fasta
# 2) list of positions for which to calculate complexity in bed format
# 3) length of the regions for which to calculate complexity (eg 10 calculates complexity for a region + and - 10bp from a mutation


#span = '10'
#fasta = 'sample_files/spades07_BM03_scaffolds_sample.fasta'
#bed = 'sample_files/intervals.bed'

span = sys.argv[1]
fasta = sys.argv[2]
bed = sys.argv[3]



def findReps(reg,reps):
    P = {}
    for ele in reps:
        pat = re.compile('('+ele+'){2,}')
        m = pat.finditer(reg)
        for match in m:
            span = match.span()
            group = match.group()
            P[span] = group, len(group)
    
    P_span = []
    if P.keys():
        spans = sorted(list(P.keys()))
        keep = [spans[0]]
        for i in range(len(spans))[1:]:
            prev0 = spans[i-1]
            prev = keep[-1]
            ele = spans[i]
            if (set(range(ele[0],ele[1])).intersection(set(range(prev[0],prev[1])))!=set()) & (set(list(P[ele][0])) == set(list(P[prev0][0]))):
            #if set(range(ele[0],ele[1])).intersection(set(range(prev[0],prev[1]))):
                prev_len = P[prev0][1]
                ele_len = P[ele][1]
                #!!! keep only longer
                if prev_len>=ele_len:
                    #print('shorter')
                    continue
                else:
                    #print('longer')
                    keep = keep[:-1]
                    keep.append(spans[i])

            elif (set(range(ele[0],ele[1])).intersection(set(range(prev[0],prev[1])))!=set()) & (set(list(P[ele][0])) != set(list(P[prev0][0]))):
                newrang = set(range(ele[0],ele[1])) - set(range(prev[0],prev[1]))
                nstart = min(newrang)
                nstop = max(newrang)
                newspan = nstart, nstop
                #print('overlap but diff')
                keep.append(newspan)
            else:
                #print('no overlap')
                keep.append(spans[i])
        R = {}
        for ele in keep:
            rep_reg = reg[ele[0]:ele[1]]
            R[ele] = rep_reg
    else:
        R = {}
    return(R)



# function for finding sequence complexity in genetic sequence
def calcComplexity(reg,max_rep):
    #max_rep = 10
    #reg = region

    T = [reg]
    uniquelen = []

    for rep in list(reversed(range(1,max_rep+1))):
        reps = [''.join(ele) for ele in list(permutations(['A','C','T','G'],rep))]
        #print(reps)
        S = []
        for sequence in T:
            findrep = findReps(sequence,reps)
            uniquelen.append(len(findrep)*rep)
            # splitting region by repeats
            rep_span = [0]
            for ele in sorted(list(findrep.keys())):
                rep_span+=[ele[0]]
                rep_span+=[ele[1]]
            rep_span.append(len(sequence))
            i = 0
            for j in range(int(len(rep_span)/2)):
                seq = sequence[rep_span[i]:rep_span[i+1]]
                i+=2
                if seq:
                    S.append(seq)
        T=S
        #print(T)

    uni_length = sum(uniquelen) + sum([len(ele) for ele in T])
    tot_length = len(region)
    return(uni_length/tot_length)


def readFasta(fpath):
    ref = SeqIO.parse(fpath,'fasta')
    f = {record.id:str(record.seq) for record in ref}
    return(f)


def readMuts(bedpath):
    fh = open(bedpath,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    return(k)


if __name__ == '__main__':
    
    hits = readMuts(bed)
    genome = readFasta(fasta)
    span_length = int(span)
    T = []
    Reg = []
    C = []
    wh = open('calcComplexity.out','w')
    wh.write('chrom\tstart\tstop\tcomplexity\n')
    for hit in hits:
        chrom = hit[0]
        pos0 = int(hit[1])
        pos1 = int(hit[2])
        ref_chrom = genome[chrom]
        tri = ref_chrom[pos0-1:pos0+2]
        region = ref_chrom[pos0-1-span_length:pos0+2+span_length]
        T.append(tri)
        Reg.append(region)
        complexity = calcComplexity(region, span_length)
        f = '\t'.join([hit[0],hit[1],hit[2],str(complexity)])
        wh.write(f+'\n')
    wh.flush()
    wh.close()

"""calculateDxy.py: A scikit-allel wrapper to calculate dxy in windows and per chromosome/contig.
   To get an unbiased dxy, generate vcf file including non-variant sites."""

__author__ = "anna fijarczyk"


import allel; print('Using scikit-allel', allel.__version__)
import pandas as pd
import numpy as np
from collections import defaultdict
import sys, getopt


#input_vcf = 'sample.vcf.gz'
#input_pop1 = 'pop1_mosaicbeer.tab'
#input_pop2 = 'pop2_4eflv2.tab'
#winsize = 10000
#stepsize = 5000

printhelp = 'calculateDxy.py -i <vcf> -p <tab> -r <tab> -w <int> -s <int>\n\
                \t-i [--input]: <file> vcf file (indexed with tabix, can be bgzipped) [required]\n\
                \t-p [--population1]: <file> 1-column file with sample names from pop 1 [required]\n\
                \t-r [--population2]: <file> 1-column file with sample names from pop 2 [required]\n\
                \t-w [--window-size]: <int> size of the window [10000]\n\
                \t-s [--step-size]: <int> step size [5000]\n'

def main(argv):
    input_vcf = ''
    input_pop1 = ''
    input_pop2 = ''
    winsize = ''
    stepsize = ''
    
    try:
        opts, args = getopt.getopt(argv,"i:p:r:w:s:",["input=","population1=","population2=","window-size=","step-size="])
    except getopt.GetoptError:
        print(printhelp)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(printhelp)
            sys.exit()
        elif opt in ("-i", "--input"):
            input_vcf = arg
        elif opt in ("-p", "--population1"):
            input_pop1 = arg
        elif opt in ("-r", "--population2"):
            input_pop2 = arg
        elif opt in ("-w", "--window-size"):
            winsize = arg
        elif opt in ("-s", "--step-size"):
            stepsize = arg
    opt_list = [ele[0] for ele in opts]
    if ("-w" not in opt_list) and ("--window-size" not in opt_list):
        winsize = "10000"
    if ("-s" not in opt_list) and ("--step-size" not in opt_list):
        minlength = "5000"
    if ("-i" not in opt_list) and ("--input" not in opt_list):
        print("Missing input vcf")
        print(printhelp)
        sys.exit(2)
    if ("-p" not in opt_list) and ("--population1" not in opt_list):
        print("Missing samples from population 1")
        print(printhelp)
        sys.exit(2)
    if ("-r" not in opt_list) and ("--population2" not in opt_list):
        print("Missing samples from population 2")
        print(printhelp)
        sys.exit(2)

    return(input_vcf, input_pop1, input_pop2, winsize, stepsize)
    #return opts


# Function for generating windows array [start, stop] for the chromosome (window will end with the last variant)
def getWindows(varlist, window, step):
    start = 1
    end = varlist[-1]
    W = []
    while (start+step) < end:
        stop = start + window
        W.append([start,stop])
        start+=step
    W.append([start,end])
    aW = np.array(W)
    return(aW)

def getDxy(vcf,pop1,pop2,window_size,step_size):

    # Getting the samples
    fh1 = open(pop1,'r').readlines()
    spop1 = [(ele.split()[0],'pop1') for ele in fh1]
    fh2 = open(pop2,'r').readlines()
    spop2 = [(ele.split()[0],'pop2') for ele in fh2]
    pops = spop1 + spop2
    Pops = {a:b for a,b in pops}
    Samples = list(Pops.keys())
    #Samples

    print("Reading vcf")
    callset = allel.read_vcf(vcf,['samples','variants/CHROM','variants/POS','calldata/GT'],samples=Samples)
    samples = callset['samples']
    chromosomes = callset['variants/CHROM']
    positions = callset['variants/POS']
    gts = callset['calldata/GT']
    variants = callset['variants/POS']
    idx = allel.ChromPosIndex(chromosomes, positions)
    chroms = []
    for cr in chromosomes:
        if cr not in chroms:
            chroms.append(cr)
    #chroms = list(set(chromosomes))

    # Getting sample indices
    populations = []
    for ele in samples:
        if ele in Pops.keys():
            populations.append(Pops[ele])
        else:
            populations.append('other')
    ds = pd.DataFrame({'sample':samples,'pop':populations})
    samples_callset_index = [list(samples).index(s) for s in ds['sample']]
    ds['callset_index'] = samples_callset_index
    dpops = defaultdict(list)
    for a,b in ds[['pop','callset_index']].values.tolist():
        dpops[a].append(b)

    print("Calculating per-window dxy")
    Dxy = []
    Glob = []
    for chrom in chroms:
        #print(chrom)
        chr_slice = idx.locate_key(chrom)
        # Getting genotypes
        chr_gts = gts[chr_slice]
        # Getting windows
        chr_vars = variants[chr_slice]
        win_array = getWindows(chr_vars,int(window_size), int(step_size))
        # Getting allele counts
        g = allel.GenotypeArray(chr_gts)
        counts_pop1 = g.count_alleles_subpops(dpops,max_allele=1)['pop1']
        counts_pop2 = g.count_alleles_subpops(dpops,max_allele=1)['pop2']
        # Calculating dxy in windows
        dxy, windows, nbases, counts = allel.windowed_divergence(chr_vars,counts_pop1,counts_pop2,windows=win_array)
        # Putting results in a dataframe
        chrDxy = pd.DataFrame({'dxy':dxy,'start':windows[:,0],'stop':windows[:,1],'counts':counts})
        chrDxy['mid'] = chrDxy['start'] + (chrDxy['stop'] - chrDxy['start'])/2
        chrDxy['chrom'] = chrom
        Dxy.append(chrDxy)
        # Calculating global dxy
        global_dxy = allel.sequence_divergence(chr_vars,counts_pop1,counts_pop2)
        Glob.append(global_dxy)
    dDxy = pd.concat(Dxy).reset_index().drop(columns=['index'])
    dDxy['window'] = list(range(len(dDxy['dxy'])))
    dDxy.to_csv('calculateDxy.out',sep='\t',index=False,header=True,columns=['chrom','start','stop','mid','window','counts','dxy'])
    #return(dDxy)
    dG = pd.DataFrame({'chrom':chroms,'dxy':Glob})
    dG.to_csv('calculateDxy_global.out',sep='\t',index=False,header=True)


if __name__ == '__main__':
    input_vcf, input_pop1, input_pop2, winsize, stepsize = main(sys.argv[1:])
    getDxy(input_vcf, input_pop1, input_pop2, winsize, stepsize)
    
    

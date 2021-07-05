"""calculatePCADist.py: A scikit-allel wrapper to calculate distance of a sample/samples from the PCA centroid 
   relative to other group of samples in non-overlapping windows of predetermined number of SNPs.
   Only sites with segreating biallelic variants and no missing genotypes are taken into account."""

__author__ = "anna fijarczyk"


import allel; print('Using scikit-allel', allel.__version__)
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.spatial import distance
import sys, getopt


#input_vcf = 'sample.vcf.gz'
#input_pop1 = 'pop1_mosaicbeer.tab'
#input_pop2 = 'pop2_4eflv2.tab'
#winsize = 10000

printhelp = 'calculatePCADist.py -i <vcf> -p <tab> -r <tab> -w <int>\n\
                \t-i [--input]: <file> vcf file (indexed with tabix, can be bgzipped) [required]\n\
                \t-p [--population1]: <file> 1-column file with sample names from the background population [required]\n\
                \t-r [--population2]: <file> 1-column file with sample name(s) from the focal population [required]\n\
                \t-w [--window-size]: <int> size of the window (number of SNPs) [100]\n'

def main(argv):
    input_vcf = ''
    input_pop1 = ''
    input_pop2 = ''
    winsize = ''
    
    try:
        opts, args = getopt.getopt(argv,"i:p:r:w:",["input=","population1=","population2=","window-size="])
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
    opt_list = [ele[0] for ele in opts]
    if ("-w" not in opt_list) and ("--window-size" not in opt_list):
        winsize = "100"
    if ("-i" not in opt_list) and ("--input" not in opt_list):
        print("Missing input vcf")
        print(printhelp)
        sys.exit(2)
    if ("-p" not in opt_list) and ("--population1" not in opt_list):
        print("Missing samples from the background population")
        print(printhelp)
        sys.exit(2)
    if ("-r" not in opt_list) and ("--population2" not in opt_list):
        print("Missing sample(s) from the focal population")
        print(printhelp)
        sys.exit(2)

    return(input_vcf, input_pop1, input_pop2, winsize)
    #return opts


# Function for calculating distance of the sample from the centroid
def runPCA(genotypes,**kwargs):
    pca = allel.pca(genotypes,n_components=2)[0]
    centroid = np.mean(pca,axis=0)
    df1 = pd.DataFrame(pca,columns=['x','y'])
    #dst = distance.euclidean(pca[0], centroid)
    df1['dst'] = df1.apply(lambda x: distance.euclidean([x['x'],x['y']],centroid),axis=1)
    indices_pop1 = kwargs['pop_1']
    indices_pop2 = kwargs['pop_2']
    group_mean = np.mean(df1.loc[indices_pop1,'dst'])
    focal_mean = df1.loc[indices_pop2,'dst']
    dist = focal_mean/group_mean
    return(dist)

def getPCADist(vcf,fpop1,fpop2,window_size):
    
    # Getting the samples
    fh1 = open(fpop1,'r').readlines()
    spop1 = [(ele.split()[0],'pop1') for ele in fh1]
    fh2 = open(fpop2,'r').readlines()
    spop2 = [(ele.split()[0],'pop2') for ele in fh2]
    pops = spop1 + spop2
    Pops = {a:b for a,b in pops}
    Samples = list(Pops.keys())

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

    print("Calculating pop distance from the centroid")
    Dist = []
    for chrom in chroms:
        #print(chrom)

        chr_slice = idx.locate_key(chrom)
        chr_vars = variants[chr_slice]
        # Getting genotypes
        chr_gts = gts[chr_slice]
        chr_gts

        # Filtering out rows (positions) with missing genotypes
        missing = allel.GenotypeArray(chr_gts).is_missing()
        bool_missing = missing.any(axis=1)
        chr_nomissing = chr_gts[~bool_missing]
        chr_nomissing
        chr_vars_nomissing = chr_vars[~bool_missing]

        # Retaining rows (positions) with segregating genotypes
        segs = allel.GenotypeArray(chr_nomissing).count_alleles()>0
        bool_segs = segs.all(axis=1)
        chr_segregating = chr_nomissing[bool_segs]
        chr_vars_segregating = chr_vars_nomissing[bool_segs]
        #chr_segregating.shape, chr_vars_segregating.shape

        # Converting genotypes to one code number
        #chr_nalt = allel.GenotypeArray(chr_gts).to_n_alt(fill=-1)
        chr_nalt = allel.GenotypeArray(chr_segregating).to_n_alt()
        chr_nalt.shape

        ### This is optional - locating unlinked variants
        #unlink = allel.locate_unlinked(chr_nalt, size=100, step=50, threshold = 0.1)
        #chr_unlink = chr_nalt[unlink]
        #chr_vars_unlink = chr_vars_segregating[unlink]

        # Calculating distance
        win_stat = allel.moving_statistic(chr_nalt,runPCA,size=int(window_size),pop_1=dpops['pop1'],pop_2=dpops['pop2'])
        flat_stat = np.concatenate(win_stat)
        starts = chr_vars_segregating[0:len(chr_vars_segregating):int(window_size)]
        stops = chr_vars_segregating[int(window_size)-1:len(chr_vars_segregating):int(window_size)]
        wf = pd.DataFrame({'chrom' : chrom, 'dist' : flat_stat, 'SNP_start' : starts[:len(flat_stat)], 
                            'SNP_stop' : stops[:len(flat_stat)],'SNPs' : int(window_size)})
        Dist.append(wf)

    dW = pd.concat(Dist)
    dW['mid'] = dW['SNP_start'] + (dW['SNP_stop'] - dW['SNP_start'])/2
    dW['window'] = list(range(len(dW['dist'])))
    dW.to_csv('calculatePCADist.out',sep='\t',index=False,header=True,columns=['chrom','SNP_start','SNP_stop','mid','window','SNPs','dist'])
    #return(dW)


if __name__ == '__main__':
    input_vcf, input_pop1, input_pop2, winsize = main(sys.argv[1:])
    getPCADist(input_vcf, input_pop1, input_pop2, winsize)
    
    

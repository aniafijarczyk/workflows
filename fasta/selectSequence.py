from Bio import SeqIO
import sys, getopt

#c = 'exclude'
#f = 'contigs_BM03.txt'
#l = 10000
#i = '../03_finalassemblies_raw/spades07_BM03_scaffolds.fasta'


printhelp =  'selectSequence.py -t <type> -c <contigs> -l <length> -f <fasta>\n\
                \t-t [--type]: what to do with contigs - include or exclude [exclude]\n\
                \t-c [--contigs]: text file with the list of contigs in one column, optional\n\
                \t-l [--length]: min length of contigs to include [1000]\n\
                \t-f [--fasta]: fasta file with contigs, required\n'

#print(printhelp)

def main(argv):
    itype = ''
    contigs = ''
    minlength = ''
    inputfasta = ''
    
    try:
        opts, args = getopt.getopt(argv,"t:c:l:f:",["type=","contigs=","length=","fasta="])
    except getopt.GetoptError:
        print(printhelp)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(printhelp)
            sys.exit()
        elif opt in ("-t", "--type"):
            itype = arg
        elif opt in ("-c", "--contigs"):
            contigs = arg
        elif opt in ("-l", "--length"):
            minlength = arg
        elif opt in ("-f", "--fasta"):
            inputfasta = arg
    opt_list = [ele[0] for ele in opts]
    if ("-t" not in opt_list) and ("--type" not in opt_list):
        itype = "exclude"
    if ("-l" not in opt_list) and ("--length" not in opt_list):
        minlength = 1000
    if ("-c" not in opt_list) and ("--contigs" not in opt_list):
        contigs = []
    if ("-f" not in opt_list) and ("--fasta" not in opt_list):
        print(printhelp)    
        sys.exit(2)
    print('Fasta file ',inputfasta)
    if contigs:
        print('Contigs file ',contigs)
    else:
        print('No contig file, keeping all contigs')
    print('What to do with contigs? ',itype)
    print('Keeping contings with min length',minlength)
    return(itype,contigs,minlength,inputfasta)
    #return opts



def filterSeq(fasta, contigs, ctype, minlength):
    if contigs:
        fh = open(contigs,'r').readlines()
        k = [ele.split()[0] for ele in fh]
    else:
        k = []

    fastaID = fasta.split('/')[-1].replace('.fasta','').replace('.fna','').replace('.fa','').replace('.fas','')
    #print('Name of fasta ',fastaID)

    R = []
    for record in SeqIO.parse(fasta,'fasta'):
        if ctype == 'include':
            if record.id in k:
                if len(record.seq) >= minlength:
                    R.append(record)

        elif ctype == 'exclude':
            if record.id not in k:
                if len(record.seq) >= minlength:
                    R.append(record)

    #return(R)
    SeqIO.write(R,fastaID+'_selected.fasta','fasta')


if __name__ == '__main__':
    exclude,contigs,minlength,inputfasta = main(sys.argv[1:])
    #print(exclude,contigs,minlength, inputfasta)
    filterSeq(inputfasta,contigs,exclude,int(minlength))

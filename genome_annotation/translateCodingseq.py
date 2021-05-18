from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob
from collections import defaultdict



fastas=glob.glob('*.codingseq')
frames = 'findIncomplete_frame.out'


def readMe(fname):
  fh = open(fname,'r')
  linie = fh.readlines()
  k = [ele.split() for ele in linie]
  n = [ele + [ele[0].split('/')[-1]] for ele in k]
  d = defaultdict(list)
  for i in n:
    d[i[-1]].append(i)
  return(d)

D = {}

#wh = open('check.out','w')

for fasta in fastas:
  fastaID = fasta.replace('.codingseq','')
  
  R = []
  print(fastaID)
  for record in SeqIO.parse(fasta,'fasta'):
    length = float(len(record.seq))
    #if length/3 != abs(length/3):
    wh.write(fastaID+'\t'+record.id+'\t'+str(length)+'\t'+str(length/3)+'\n')

    if fastaID in D.keys():
      gene_frame = {g[1]:g[2] for g in D[fastaID]}
      if record.id.split('.')[0] in gene_frame.keys():
        start = int(gene_frame[record.id.split('.')[0]])
        aaseq = record.seq[start:].translate()
        print(fastaID,start,record.id.split('.')[0])
      else:
        aaseq = record.seq.translate()
    else:
      aaseq = record.seq.translate()

    newrecord = SeqRecord(aaseq,record.id,'','')
    R.append(newrecord)
  SeqIO.write(R,fastaID+'.aa','fasta')


#wh.flush()
#wh.close()

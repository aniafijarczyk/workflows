import glob
from collections import defaultdict
import re


gff=glob.glob('*.gff')

wh = open('findIncomplete.out','w')
wh2 = open('findIncomplete_frame.out','w')

for gff in gffs:
  fh = open(gff,'r')
  linie = fh.readlines()
  k = [ele.split() for ele in linie if ele[0] not in '#']
  # excluding 'gene' & selecting only paradoxus (cerevisiae doesnt have start and stop codons
  k1 = [ele for ele in k if ele[2] in ['CDS','exon','intron','mRNA','start_codon','stop_codon','transcript'] if ele[0].startswith('utg')]

  k2 = []
  for row in k1:
    dicID = {i.split('=')[0]:i.split('=')[1] for i in row[-1].split(';')}
    if row[2] in ['mRNA','transcript']:
      newid = dicID['ID']
    elif row[2] in ['exon','intron','start_codon','stop_codon','CDS']:
      newid = dicID['Parent']

    newrow = row + [newid]
    k2.append(newrow)

  # grouping sequences 
  d = defaultdict(list)
  for ele in k2:
    d[ele[-1]].append(ele)


  for ele in d.keys():
    types = [col[2] for col in d[ele]]
    if ('start_codon' not in types) or ('stop_codon' not in types):
      #print(gff,ele)
      wh.write('# '+gff+'\n')
      for row in d[ele]:
        f = '\t'.join(row)+'\n'
        wh.write(f)

    if ('start_codon' not in types):
      strand = d[ele][0][6]
      cdses = [i for i in d[ele] if i[2] == 'CDS']
      #print(ele)
      if strand=='+':
        #first_cds_frame = cdses[0][7]
        first_cds_frame = '0'
        frame = first_cds_frame
      elif strand=='-':
        #last_cds_frame = cdses[-1][7]
        last_cds_frame = '0'
        frame = last_cds_frame
      wh2.write(gff+'\t'+ele+'\t'+frame+'\n')    


wh.flush()
wh.close() 

wh2.flush()
wh2.close()
  

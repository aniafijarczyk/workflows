import glob

insert_size = 280
seed_input = "../../../DATA/assemblies/CM001753.1_ophniu10.fasta"
samples = ['O_triangulosporum','HP32','O_populinum','O_quercus','O_montium_2340']

#flist = glob.glob("/project/biosafe/DATA/DED/world-sequencing/seqcleaned/HI.4257.001.BioOHT_*R1P*")
#flist = glob.glob("../../../DATA/reads/varspecies_trimmomatic/HI*_outR1P.fastq.gz")
#IDs = [(ele.split("/")[-1].split('.')[4].replace("_outR1P",""),ele) for ele in flist if ]
#flistDict = {a:b for a,b in zip(IDs,flist)}

flistDict = {ele : glob.glob('../../../DATA/reads/varspecies_trimmomatic/HI*'+ele+'_outR1P.fastq.gz')[0] for ele in samples}
print(flistDict)

for ids in flistDict.keys():
  wh = open("config_"+ids+".txt","w")
  wh.write("Project name         = "+ids+"_mit\n")
  wh.write("Insert size          = "+str(insert_size)+"\n")
  wh.write("Insert size aut      = yes\n")
  wh.write("Read Length          = 150\n")
  wh.write("Type                 = mito\n")
  wh.write("Genome Range         = 10000-90000\n")
  wh.write("K-mer                = 39\n")
  wh.write("Insert Range         = 1.6\n")
  wh.write("Insert Range strict  = 1.2\n")
  wh.write("Single/Paired        = PE\n")
  wh.write("Max memory           = 30\n")
  wh.write("Extended log         = 0\n")
  wh.write("Save assembled reads = yes\n")
  wh.write("Combined reads       =\n")
  wh.write("Forward reads        = "+flistDict[ids]+"\n")
  wh.write("Reverse reads        = "+flistDict[ids].replace("_outR1P","_outR2P")+"\n")
  wh.write("Seed Input           = "+seed_input+"\n")
  wh.write("Reference            =\n")
  wh.write("Chloroplast sequence =\n")
  
  wh.flush()
  wh.close()


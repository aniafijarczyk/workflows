### Genome assembly from illumina PE reads for medium sized genomes (fungi)

[1. Quality check (fastqc)](#1-Quality-check)  
[2. Adapter trimming (trimmomatic)](#2-Adapter-trimming)  
[3. *(optionally)* Read merging (bbmerge)](#3-Read-merging)  
[4. Assembly (spades,abyss)](#4-Assembly)  
[5. Assembly evaluation (quast)](#5-Assembly-evaluation)  
[6. Contamination](#6-contamination)  
[7. Mitochondrial genome assembly](#7-mitochondrial-genome-assembly)  
[8. Contig coverage](#8-contig-coverage)  
[9. Gene completeness (busco)](#9-gene-completeness)  


##### 1 Quality check

Checking quality of fastq reads with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.11.8)

```
fastq=*.fastq
list=$(ls $fastq)
fastqc -o ./ ${list}
```
Check [help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for interpretation

##### 2 Adapter trimming

Using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) (v0.36)

```
for sample in AT146 Yu16
  do
  left=reads/${sample}*_R1.fastq.gz
  right=reads/${sample}*_R2.fastq.gz
  java -jar -Xmx4g /prg/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
  reads/${sample}_R1.fastq.gz reads/${sample}_R2.fastq.gz \
  trimmed/${sample}_outR1P.fastq trimmed/${sample}_outR1U.fastq trimmed/${sample}_outR2P.fastq trimmed/${sample}_outR2U.fastq \
  ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 MINLEN:21
  done
```

or with parallel
```
parallel -j 10 java -jar -Xmx4g /prg/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 reads/{}_R1.fastq.gz reads/{}_R2.fastq.gz \
trimmed/{}_outR1P.fastq trimmed/{}_outR1U.fastq trimmed/{}_outR2P.fastq trimmed/{}_outR2U.fastq \
ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 MINLEN:21 ::: $(ls -1 reads/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | cut -d'/' -f2)
```

Location of adapter sequences: /prg/trimmomatic/0.36/adapters/

##### 3 Read merging

Merging overlapping reads with bbmerge from [bbmap](https://github.com/BioInfoTools/BBMap)
```
bbmerge.sh -t=8 -in1=reads_outR1P.fastq -in2=reads_outR2P.fastq -out=reads_merged.fastq -outu1=reads_R1um.fastq -outu2=reads_R2um.fastq
```

##### 4 Assembly

Assembly of PE reads using [SPAdes](https://github.com/ablab/spades) with memory limit of 100G and 8 threads (v3.9.1)
```
spades.py -k 21,33,55,77,99 --careful --pe1-1 reads_R1P.fastq --pe1-2 reads_R2P.fastq -o spades_v1 -t 8 -m 100
```
Assembly of PE and merged reads with SPAdes
```
spades.py -k 21,33,55,77,99 --careful --pe1-1 reads_R1um.fastq --pe1-2 reads_R2um.fastq --s1 reads_merged.fastq -o spades_v2 -t 8 -m 100
```
Assembly of PE reads with [ABySS](https://github.com/bcgsc/abyss) (v2.0.2) with min base quality 15 and 64 kmer
```
R1=reads_R1P.fastq
R2=reads_R2P.fastq
abyss-pe -C abyss_v1 name=tetropii k=64 q=15 in="$R1 $R2"
```
Assembly of PE and merged reads with ABySS
```
R1=reads_R1um.fastq
R2=reads_R2um.fastq
RS=reads_merged.fastq
abyss-pe -C abyss_v2 name=tetropii k=64 q=15 in="$R1 $R2" se="$RS"
```
##### 5 Assembly evaluation
Using [quast](https://github.com/ablab/quast) (v4.3)
```
quast.py -o quast_spades -m 200 -l spades_v1,spades_v2,abyss_v1 ./spades_v1/scaffolds.fasta ./spades_v2/scaffolds.fasta ./abyss_v1/tetropii.fasta
quast.py -o quast_spades_broken -m 200 -l spades_v1,spades_v2,abyss_v1 ./spades_v1/scaffolds.fasta ./spades_v2/scaffolds.fasta ./abyss_v1/tetropii.fasta
```
##### 6 Contamination
A simple approach to filter out contaminated contigs. In short, in case of assmbling fungi, contigs are blasted against fungal proteins from SwissProt and separately against bacterial and viral proteins. Contigs with higher overlap of bac/vir than fungal proteins are marked as potential contaminants.  
Blasting contigs using [diamond](https://github.com/bbuchfink/diamond) (v0.9.29)
```
dbdir=uniprot_sprot.fasta
diamond blastx -p 10 -d $dbdir -q scaffolds.fasta -o blastx_sample_uniprot_fungi.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles staxids \
-e 0.001 --max-target-seqs 1000 --taxonlist 4751
diamond blastx -p 10 -d $dbdir -q scaffolds.fasta -o blastx_sample_uniprot_bacvir.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles staxids \
-e 0.001 --max-target-seqs 1000 --taxonlist 2,2157,10239
```
Find overlapping hits in each output
```
getBestOverlapFromBlastx.py fungi blastx_sample_uniprot_fungi.tab
getBestOverlapFromBlastx.py bacvir blastx_sample_uniprot_bacvir.tab
```
Combine files with fungal and bacvir overlaping hits
```
cat blastx_sample_uniprot_fungi_bestoverlap.tab blastx_sample_uniprot_bacvir_bestoverlap.tab > blastx_sample_merged.tab
```
Find best hit from fungal and bacvir hits
```
getBestOverlapFromMerged.py blastx_sample_merged.tab
```
Compare fungal and bacvir hits and find contaminated contigs
```
python calculateOverlap.py blastx_sample_merged_bestoverlap.tab spades_TB5_clean_softmasked.fasta.fai
``` 

##### 7 Mitochondrial genome assembly

Preparing config files (seed mtDNA genome, insert length, path to fastq files) 
```
python makeConfigFile.py
```
Assembling mtDNA with [Novoplasty](https://github.com/ndierckx/NOVOPlasty) (v3.8.3)
```
configFiles=config_*.txt
for config in $configFiles
  do
  perl $DIR/NOVOPlasty3.8.3.pl -c $config
  done

```
Aligning mtDNA to whole assembly to identify mtDNA contigs with [nucmer](http://mummer.sourceforge.net/) (mummer v4.0)
```
nucmer --maxgap=1000 --mincluster=100 --prefix=nucmer_mtDNA scaffolds.fasta mtdna_assembly.fasta
delta-filter -q -r nucmer_mtDNA.delta > nucmer_mtDNA_qr.filter
show-coords -rclT nucmer_mtDNA_qr.filter > nucmer_mtDNA_qr.coords
```
Calculate % contig length covered by mtDNA
```
python getMtdnaInfo.py
```

##### 8 Contig coverage

Checking read depth of contigs to exclude low-covered contigs (eg ncov < 0.05) 
Mapping reads to assembly and calculating depth
```
ref=scaffolds.fasta
reads1=reads_outR1P.fastq.gz
reads2=reads_outR2P.fastq.gz
bwa mem -t 3 ${ref} ${reads1} ${reads2} | samtools sort -o bams/scaffolds.bam -
samtools index bams/scaffolds.bam
samtools depth -a bams/scaffolds.bam > depths/scaffolds_depth.out
```
Calculating mean read depth per contig (cov_mean) and normalized read depth (ncov = cov_mean/median assembly)
```
python getCoverage.py
```

##### 9 Gene completeness
Finding gene completeness with [busco](https://busco.ezlab.org/) (v3.0.1)
```
lin=sordariomyceta_odb9
run_BUSCO.py -i scaffolds_${sample}.fasta -m geno -o busco_${sample} -l $lin -sp magnaporthe_grisea -c 4
or
run_BUSCO.py -i scaffolds_${sample}.fasta -m geno -o busco_${sample} -l $lin -sp magnaporthe_grisea -c 1
```
Summarize info
```
python getBusco.py
```

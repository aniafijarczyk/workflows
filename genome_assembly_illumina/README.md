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

Checking quality of fastq reads with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
fastq=*.fastq
list=$(ls $fastq)
fastqc -o ./ ${list}
```
Check [help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) for interpretation

##### 2 Adapter trimming

Using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

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


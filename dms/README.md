### Analysis of DMS data


[1. Quality check of fastq](#1-quality-check-of-fastq)  
[2. Read trimming](#2-read-trimming)  
[3. Mapping](#3-mapping)  
[4. Analyse counts](#4-analyse-counts)



##### 1 Quality check of fastq

Using [fastqc/0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```
fastq=./DATA/*.fastq
list=$(ls $fastq)
fastqc -o ./ ${list}
```
[MultiQC](https://multiqc.info/)



##### 2 Read trimming

Trimming adapters & low quality bases with [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) in a loop
```
for sample in AT146 Yu16
  do
  left=reads/${sample}*_R1.fastq.gz
  right=reads/${sample}*_R2.fastq.gz
  java -jar -Xmx4g /prg/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
  reads/${sample}_R1.fastq.gz reads/${sample}_R2.fastq.gz \
  trimmed/${sample}_outR1P.fastq trimmed/${sample}_outR1U.fastq trimmed/${sample}_outR2P.fastq trimmed/${sample}_outR2U.fastq \
  ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:3:24 MINLEN:21
  done
```

or with parallel
```
parallel -j 10 java -jar -Xmx4g /prg/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 reads/{}_R1.fastq.gz reads/{}_R2.fastq.gz \
trimmed/{}_outR1P.fastq trimmed/{}_outR1U.fastq trimmed/{}_outR2P.fastq trimmed/{}_outR2U.fastq \
ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:3:24 MINLEN:21 ::: $(ls -1 reads/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | cut -d'/' -f2)
```
Location of adapter sequences: /prg/trimmomatic/0.36/adapters/


Trimming adapters and low quality bases with [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html)



##### 3 Mapping

Mapping reads locally with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)



##### 4 Analyse counts

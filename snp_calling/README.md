### Scheme for SNP calling from aligned short read illumina reads


1. Quality check of fastq (fastqc, multiqc)
2. Adapter trimming (trimmomatic, cutadapt)
3. Mapping (bwa mem)
4. *(optional)* Checking mapping quality (goleft)
5. *(optional)* Generating read depth per site (samtools)
6. Removing duplicates (picard-tools)
7. SNP calling with samtools/bcftools
8. SNP calling with GATK Haplotype Caller
9. SNP calling with Freebayes



##### 1. Quality check of fastq 

```
fastq=../../../DATA/rnaseq/*.fastq
list=$(ls $fastq)
fastqc -o ./ ${list}
```

##### 2. Adapter trimming


```
for sample in AT146 Yu16
  do
  left=${sample}*_R1.fastq.gz
  right=${sample}*_R2.fastq.gz
  java -jar -Xmx4g /prg/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33
  ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz 
  output/${sample}_outR1P.fastq output/${sample}_outR1U.fastq output/${sample}_outR2P.fastq output/${sample}_outR2U.fastq 
  ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 MINLEN:21
  done
```

```
parallel -j 10 java -jar -Xmx4g /prg/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 {}_R1.fastq.gz {}_R2.fastq.gz 
output/{}_outR1P.fastq output/{}_outR1U.fastq output/{}_outR2P.fastq output/{}_outR2U.fastq 
ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 MINLEN:21 ::: $(ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//')
```

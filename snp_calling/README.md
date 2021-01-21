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


In a loop
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

With parallel
```
parallel -j 10 java -jar -Xmx4g /prg/trimmomatic/0.33/trimmomatic-0.33.jar PE -phred33 reads/{}_R1.fastq.gz reads/{}_R2.fastq.gz \
trimmed/{}_outR1P.fastq trimmed/{}_outR1U.fastq trimmed/{}_outR2P.fastq trimmed/{}_outR2U.fastq \
ILLUMINACLIP:TruSeq3-PE-adapters.fa:6:20:10 MINLEN:21 ::: $(ls -1 reads/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | cut -d'/' -f2)
```

Location of adapter sequences: /prg/trimmomatic/0.36/adapters/


##### 3. Mapping

Mapping & sorting, indexing and calculating stats with 2 threads (bwa/0.7.17, samtools/1.8)
```
bwa index reference.fasta
samtools faidx reference.fasta

for sample in AT146 Yu16
  do
  if [ ! -f bams/"${sample}_SpA_sorted.bam" ]
  then
    bwa mem -t 2 reference.fasta ${sample}_outR1P.fastq ${sample}_outR2P.fastq | samtools sort -@ 2 -o bams/${sample}_SpA_sorted.bam -
    samtools index bams/${sample}_SpA_sorted.bam
    samtools stats bams/${sample}_SpA_sorted.bam > bams/${sample}_SpA_sorted_stat
  fi
  done

```


##### 4. Checking mapping quality

Mapping coverage from bam index files with [goleft](https://github.com/brentp/goleft)
```
goleft indexcov --directory stats_indexcov/ bams/*.bam
```

##### 5. Generating read depth per site

```
samtools depth -a bams/${sample}_SpA_sorted.bam > bams/${sample}_SpA_sorted_depth.out
```

##### 6. Removing duplicate reads


Removing duplicates with picard (picard-2.18)
```
java -jar ${dir}/picard/build/libs/picard.jar MarkDuplicates I=${sample}_SpA_sorted.bam O=${sample}_SpA_rmdup.bam M=${sample}_SpA_picard.metrics REMOVE_DUPLICATES=true
```

Creating file with commands
```
rm queue.txt
touch queue.txt
for sample in bams/*.bam
  do
  nameis=$(echo $sample | sed 's/_sorted.bam//' | cut -d"/" -f 2)
  echo $nameis
  echo 'java -jar '${dir}'/picard/build/libs/picard.jar MarkDuplicates I=bams/'${nameis}'_sorted.bam O='${nameis}'_rmdup.bam M='${nameis}'_picard.metrics REMOVE_DUPLICATES=true' >> queue.txt
  done
```

##### 7. SNP calling with samtools/bcftools

Realigning reads around indels
```
for bam in bams/*.bams
  do
  nameis=$(echo $bam | sed 's/.bam//g')
  echo $nameis
  samtools calmd -bAr -@2 ${bam} reference.fasta > ${nameis}_calmd.bam
  done
```

SNP calling with [bcftools mpileup](http://samtools.github.io/bcftools/bcftools.html#mpileup) and [bcftools call](http://samtools.github.io/bcftools/bcftools.html#call)
[I think bcftools mpileup doesn't work with bcftools call -m (using -c instead)]
```
ls bams/*_calmd.bam > bam_list
bcftools mpileup -C50 -f reference.fasta -q4 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP -Ou -b bams_list | \
bcftools call -cv -f gq -Oz -o snp_bcftools.vcf.gz -
tabix -p vcf snp_bcftools.vcf.gz
```

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
  nameis=$(echo $sample | sed 's/_sorted.bam//' | cut -d'/' -f 2)
  echo $nameis
  echo 'java -jar '${dir}'/picard/build/libs/picard.jar MarkDuplicates I=bams/'${nameis}'_sorted.bam O='${nameis}'_rmdup.bam M='${nameis}'_picard.metrics REMOVE_DUPLICATES=true' >> queue.txt
  done
```

##### 7. SNP calling with samtools/bcftools

Realigning reads around indels with [samtools calmd](http://www.htslib.org/doc/samtools-calmd.html)
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

Some filtering with [vcf-annotate](https://manpages.debian.org/testing/vcftools/vcf-annotate.1.en.html)

```
zcat snp_bcftools.vcf.gz | grep -v "^#" | awk '{print $1"\t"$2"\t"$2"\tsnp_"NR}' > snp_bcftools_annotations
bgzip snp_bcftools_annotations
tabix -s 1 -b 2 -e 3 snp_bcftools_annotations.gz
```
```
export PERL5LIB=/home/anna/github/vcftools/src/perl
bcftools view snp_bcftools.vcf.gz | vcf-annotate -f +/d=10/D=10000/q=20/Q=10/-W/w=10 \
--fill-HWE --fill-type -n -a snp_bcftools_annotations.gz -c CHROM,FROM,TO,INFO/SNP_ID \
-d key=INFO,ID=SNP_ID,Number=1,Type=Integer,Description='SnpList' | bgzip -c > snp_bcftools_annotated.vcf.gz
```
```
StrandBias FLOAT Min P-value for strand bias (given PV4) [0.0001]
BaseQualBias FLOAT Min P-value for baseQ bias [1e-100]
MapQualBias FLOAT Min P-value for mapQ bias [0]
EndDistBias FLOAT Min P-value for end distance bias [0.0001]
a = MinAB INT Minimum number of alternate bases [2]
c = SnpCluster INT1,INT2 Filters clusters of 'INT1' or more SNPs within a run of 'INT2' bases []
D = MaxDP INT Maximum read depth [10000000]
d = MinDP INT Minimum read depth [2]
q = MinMQ INT Minimum RMS mapping quality for SNPs [10]
Q = Qual INT Minimum value of the QUAL field [10]
W = GapWin INT Window size for filtering adjacent gaps [10]
w = SnpGap INT SNP within INT bp around a gap to be filtered [10]
```

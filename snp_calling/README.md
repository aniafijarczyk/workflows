### Scheme for SNP calling from aligned short illumina reads


[1. Quality check of fastq (fastqc, multiqc)](#new-header)
2. Adapter trimming (trimmomatic, cutadapt)
3. Mapping (bwa mem)
4. *(optional)* Checking mapping quality (goleft)
5. *(optional)* Generating read depth per site (samtools)
6. Removing duplicates (picard-tools)
7. SNP calling with samtools/bcftools
8. SNP calling with GATK Haplotype Caller
9. SNP calling with Freebayes



# New Header
##### 1. Quality check of fastq 

Using (fastqc/0.11.8)
```
fastq=../../../DATA/rnaseq/*.fastq
list=$(ls $fastq)
fastqc -o ./ ${list}
```
[MultiQC](https://multiqc.info/)

##### 2. Adapter trimming

Using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) in a loop
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


##### 3. Mapping

Mapping & sorting, indexing and calculating stats with 2 threads with [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) (bwa/0.7.17)
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

##### 6. Marking duplicate reads


Marking duplicates with [picard](https://broadinstitute.github.io/picard/) (picard-2.18)
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
  echo 'java -jar '${dir}'/picard/build/libs/picard.jar MarkDuplicates I=bams/'${nameis}'_sorted.bam O=bams/'${nameis}'_rmdup.bam M=bams/'${nameis}'_picard.metrics REMOVE_DUPLICATES=true' >> queue.txt
  done
```

##### 7. SNP calling with samtools/bcftools

Realigning reads around indels with [samtools calmd](http://www.htslib.org/doc/samtools-calmd.html)
```
for bam in bams/*_rmdup.bams
  do
  nameis=$(echo $bam | sed 's/.bam//g')
  echo $nameis
  samtools calmd -bAr -@2 ${bam} reference.fasta > ${nameis}_calmd.bam
  done
```

SNP calling with [bcftools mpileup](http://samtools.github.io/bcftools/bcftools.html#mpileup) and [bcftools call](http://samtools.github.io/bcftools/bcftools.html#call) (bcftools v1.9)
[I think bcftools mpileup doesn't work with bcftools call -m (using -c instead)]
```
ls bams/*_calmd.bam > bam_list
bcftools mpileup -C50 -f reference.fasta -q4 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP -Ou -b bams_list | \
bcftools call -cv -f gq -Oz -o snp_bcftools.vcf.gz -
tabix -p vcf snp_bcftools.vcf.gz
```

Annotation of poor SNPs with [vcf-annotate](https://manpages.debian.org/testing/vcftools/vcf-annotate.1.en.html)

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
Filtering with [bcftools view](http://samtools.github.io/bcftools/bcftools.html#view) - only SNPs, with coverage min 10 per sample, and variants which pass above filters
```
bcftools view -i 'TYPE=="snp" & (DP4[0]+DP4[1]+DP4[2]+DP4[3])>=10 & FILTER=="PASS"' \
snp_bcftools_annotated.vcf.gz -Oz -o snp_bcftools.f.vcf.gz
```
Filtering with [vcftools](https://vcftools.github.io/man_latest.html) - same as above + removing low quality genotypes (GQ<20) and SNPs with >10% missing data with
```
vcftools --gzvcf snp_bcftools_annotated.vcf.gz --remove-indels --remove-filtered-all --minDP 10 --minGQ 20 --max-missing 0.9 --recode --recode-INFO-all --out snp_bcftools.f2
```

##### 8. SNP calling with GATK Haplotype Caller (gatk-4.1.0.0)

If bams generate error - readjust RG field following [this site](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)   
Here adjusting RG fields at the mapping stage (can be also done with picard - shown below under freebayes)
```
bwa mem -t 10 -R "@RG\tID:HKTFYCCXY.3\tLB:G4_1\tPL:ILLUMINA\tPU:HKTFYCCXY.3.sampleID\tSM:sampleID" ref.fasta reads/${sample}_R1P.fastq.gz reads/${sample}_R2P.fastq.gz | samtools sort -@ 10 -o bams2/${sample}_sorted.bam -
```
```
where:
# ID:flowcell.lane
# PL:ILLUMINA
# SM:sample
# LB:sample
# PU:flowcell.lane.sample
```
Getting RG info from each fastq file (may need to be adjusted to sequencing tech output): filename, sampleID, flow cell, flow cell lane, sth else
```
fastqs=reads/*_R1.fastq.gz
for fastq in $fastqs
  do
  nameis=$(echo $fastq | cut -d"/" -f2)
  core=$(echo $nameis | sed 's/_R1.fastq.gz//g')
  ids=$(echo $nameis | sed 's/_R1.fastq.gz//g' | cut -d"." -f1)
  zcat reads/${nameis} | head -n1 | awk -F":" -v name="${core}" -v var="${ids}" '{print name,var,$3,$4,$10}'
  done
```

Indexing reference
```
gatk CreateSequenceDictionary -R ref.fasta
```
Haplotype Caller
```
for bam in bams/*_rmdup.bam
  do
  sample=$(echo $bam | cut -d'/' -f2 | sed 's/_rmdup.bam//g')
  echo $sample
  gatk HaplotypeCaller \
    -R ref.fasta \
    -I $bam \
    -O gatk/${sample}.g.vcf.gz \
    -ERC GVCF \
    -G StandardAnnotation
  done
```
Combine GVCFs
```
gatk CombineGVCFs -R ref.fasta \
    -V gatk/samp1.g.vcf.gz -V gatk/samp2.g.vcf.gz -V gatk/samp3.g.vcf.gz \
    -O combined/combined.g.vcf \
    -G StandardAnnotation
done
```
Genotype GCVF
```
gatk GenotypeGVCFs -R ref.fasta \
    -V combined/combined.g.vcf \
    -O combined/combined.vcf \
    #--include-non-variant-sites true \
    -G StandardAnnotation \
    --annotation Coverage \
    --annotation QualByDepth \
    --annotation MappingQuality \
    --annotation MappingQualityRankSumTest \
    --annotation ExcessHet
```
Variant Filtration - annotating poor quality SNPs: standard hard filters + GQ < 20, DP per sample < 3  
Adjusting hard filters [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
```
gatk VariantFiltration -R ref.fasta \
    -V combined/combined.vcf \
    -O filtered/combined_VF.vcf \
    -filter "QD < 2.0" \
    --filter-name "QD" \
    -filter "MQ < 40.0" \
    --filter-name "MapQ" \
    -filter "FS > 60.0" \
    --filter-name "FS" \
    -filter "SOR > 3.0" \
    --filter-name "SOR" \
    -filter "MQRankSum < -12.5" \
    --filter-name "MQRS" \
    -filter "ReadPosRankSum < -8.0" \
    --filter-name "RPRS" \
    --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "GQ" \
    --genotype-filter-expression "DP < 3" \
    --genotype-filter-name "DP" \
    --set-filtered-genotype-to-no-call true \
    --QUIET true
```
Select Variants - only SNPs which pass above filters
```
gatk SelectVariants -R ref.fasta \
    -V filtered/combined_VF.vcf \
    -O filtered/combined_SV.vcf \
    --exclude-filtered true \
    --exclude-non-variants true \
    --remove-unused-alternates true \
    --select-type-to-exclude MIXED \
    --select-type-to-exclude SYMBOLIC \
    --select-type-to-exclude INDEL
```

##### 9. SNP calling with Freebayes (v1.3.1-17-gaa2ace8)

Readjustment of RG fields in bam files using picard, if needed
```
picard AddOrReplaceReadGroups I=bams/${sample}_rmdup.bam O=bams/${sample}_RG.bam RGID=${sampleID} RGLB=${sampleID} RGPL=ILLUMINA RGPU=${flowcell.lane}.${sampleID} RGSM=${sampleID}
```
SNP calling with [freebayes](https://github.com/freebayes/freebayes)
```
bams=bams/*_RG.bam
bamfiles_list=$(ls $bams)
freebayes -f ref.fasta -q 20 --use-best-n-alleles 4 --limit-coverage 20000 -F 0.02 -L $bamfiles_list > freebayes.vcf
```
Splitting variants (problematic in multiallelic variants)
```
$DIR/vcflib/bin/vcfallelicprimitives -kg freebayes.vcf >freebayes_calls.vcf
```
Filtering variants
```
$DIR/vcflib/bin/vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & MQM / MQMR > 0.9 & MQM / MQMR < 1.05" freebayes_calls.vcf > freebayes_calls.f.vcf
```
Selecting SNPs
```
bgzip freebayes_calls.f.vcf
tabis -p vcf freebayes_calls.f.vcf.gz
bcftools view -m2 -M2 --include 'TYPE="snp"' -O v -o freebayes_calls.f.snp.vcf freebayes_calls.f.vcf.gz
```

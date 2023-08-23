### Scheme for SNP calling from aligned short illumina reads


[1. Quality check of fastq (fastqc, multiqc)](#1-quality-check-of-fastq)  
[2. Adapter trimming (trimmomatic)](#2-adapter-trimming)  
[3. Mapping (bwa mem)](#3-mapping)  
[4. Checking mapping quality (goleft) *(optional)*](#4-checking-mapping-quality)  
[5. Checking read depth (mosdepth, samtools) *(optional)*](#5-checking-read-depth)  
[6. Marking duplicates (picard-tools)](#6-marking-duplicates)  
[7. Checking read depth again (mosdepth) *(optional)*](#7-checking-read-depth-again)  
[8. SNP calling with samtools/bcftools](#8-snp-calling-with-bcftools)  
[9. SNP calling with GATK Haplotype Caller](#9-snp-calling-with-gatk)  
[10. SNP calling with Freebayes](#10-snp-calling-with-freebayes)  



### Variables
```
fq1=ABC_1.fastq.gz
fq2=ABC_2.fastq.gz
SAMPLE=ABC
FASTQ_THREADS=8
TRIMMOMATIC_THREADS=8
PHRED=33
ADAPTERS=trimmomatic_adapters_PE.fa
MINLEN=31
REF=reference.fasta
BWA_THREADS=8
FLOWCELL_LANE=flowcell
LANE_NUMBER=lane
PLATFORM=ILLUMINA
```

### 1 Quality check of fastq

Using [fastqc v0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```
fastqc $fq1 -o ./output_${SAMPLE} -t ${FASTQC_THREADS}
```
Using [multiQC v1.14](https://multiqc.info/)
```
multiqc -f -p -o output_${SAMPLE} ./output_${SAMPLE}
```

### 2 Adapter trimming

Using [trimmomatic v0.39](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
```
trimmomatic PE -threads ${TRIMMOMATIC_THREADS} -phred${PHRED} $fq1 $fq2 \
02_trimmed/${SAMPLE}_R1P.fastq.gz 02_trimmed/${SAMPLE}_R1U.fastq.gz 02_trimmed/${SAMPLE}_R2P.fastq.gz 02_trimmed/${SAMPLE}_R2U.fastq.gz \
ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${MINLEN}
```

trimmomatic with parallel
```
parallel -j 10 java -jar -Xmx4g /prg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 reads/{}_R1.fastq.gz reads/{}_R2.fastq.gz \
02_trimmed/{}_outR1P.fastq.gz 02_trimmed/{}_outR1U.fastq.gz 02_trimmed/{}_outR2P.fastq.gz 02_trimmed/{}_outR2U.fastq.gz \
ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31 ::: $(ls -1 reads/*_1.fastq.gz | sed 's/_1.fastq.gz//' | cut -d'/' -f2)
```

Location of adapter sequences: /prg/trimmomatic/0.36/adapters/.

Script check_encoding.sh can be used to retrieve fastq phred encoding if needed.

After that goes another quality check with fastqc to compare.


### 3 Mapping

Mapping, sorting, realigning reads, indexing and calculating stats with [bwa mem v0.7.17](http://bio-bwa.sourceforge.net/bwa.shtml) and [samtools v1.17](http://www.htslib.org/doc/samtools.html)
```
bwa index ${REF}
samtools faidx ${REF}
bwa mem -t $BWA_THREADS -R "@RG\tID:"${FLOWCELL_LANE}"."${LANE_NUMBER}"\tLB:"${SAMPLE}"\tPL:"${PLATFORM}"\tPU:"${FLOWCELL_LANE}"."${SAMPLE}"\tSM:"${SAMPLE}"" \
$REF 02_trimmed/${SAMPLE}_R1P.fastq.gz 02_trimmed/${SAMPLE}_R2P.fastq.gz | \
samtools sort -@ $SAMTOOLS_THREADS -u - | \
samtools calmd -bAQr -@ $SAMTOOLS_THREADS - ${REF} > 04_bwa/${SAMPLE}_sorted.bam
samtools index 04_bwa/${SAMPLE}_sorted.bam
samtools stats 04_bwa/${SAMPLE}_sorted.bam > 04_bwa/${SAMPLE}_sorted.stat
```
RG field is readjusted following [this information](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups), where:
```
# ID:FLOWCELL_LANE.LANE_NUMBER
# PL:ILLUMINA
# SM:SAMPLE
# LB:SAMPLE
# PU:FLOWCELL_LANE.SAMPLE
```
Script check_flowcell.sh can be to retrieve flowcell lane and lane number from fastq header if needed.


### 4 Checking mapping quality (optional)

Mapping coverage from bam index files with [goleft](https://github.com/brentp/goleft)
```
goleft indexcov --directory stats_indexcov/ bams/*.bam
```

### 5 Checking read depth

Average read depth per chromosome/contig using [mosdepth v0.0.3](https://github.com/brentp/mosdepth)
```
mosdepth --fast-mode ./05_depth/${SAMPLE} ./04_bwa/${SAMPLE}_sorted.bam
```

Read depth per site with samtools (mosdepth is faster)
```
samtools depth -a 04_bwa/${SAMPLE}_sorted.bam > 04_bwa/${SAMPLE}_sorted_depth.out
```

### 6 Marking duplicates

Marking duplicate reads with [picard v2.18.29](https://broadinstitute.github.io/picard/)
```
picard MarkDuplicates I=./04_bwa/${SAMPLE}_sorted.bam O=./06_rmdup/${SAMPLE}_rmdup.bam M=./06_rmdup/${SAMPLE}_picard.metrics REMOVE_DUPLICATES=true
samtools index 06_rmdup/${SAMPLE}_rmdup.bam
```

### 7 Checking read depth again

Read depth per site and per 1000bp window with mosdepth
```
mosdepth --fast-mode --by 1000 ./07_depth/${SAMPLE} ./06_rmdup/${SAMPLE}_rmdup.bam
```


### 8 SNP calling with bcftools

Variant calling with [bcftools v1.17](https://samtools.github.io/bcftools/bcftools.html)
```
BAMS=bams_list.txt
bcftools mpileup -C50 -f ${REF} -min-MQ 4 -min-BQ 13 --skip-any-set 1796 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP -Ou -b ${BAMS} | \
bcftools call -mv -f gq -Oz -o snp_bcftools.vcf.gz -
tabix -p vcf snp_bcftools.vcf.gz

# --skip-any-set: Skip reads with any of the FLAG bits set;
# 1796: read unmapped (0x4), not primary alignment (0x100), read fails platform/vendor quality checks (0x200), read is PCR or optical duplicate (0x400)
```

Checking the stats
```
bcftools stats snp_bcftools.vcf.gz > snp_bcftools.stats
```

Example filtering:  
- removing:  
  - variants with average DP across samples < 10  
  - total sum of DP > 20000  
  - variant quality < 20  
  - read mapping quality < 40  
  - ALT allele frequency across samples < 0.05  
  - variants with allele count of 0  
- keeping:  
  - biallelic variants  
  - variants with at least 1 alternative allele  

```
bcftools filter -e "AVG(FORMAT/DP)<10 || INFO/DP>20000 || QUAL<20 || MQ<40 || (AF < 0.05) || (AC == 0)" -Ou snp_bcftools.vcf.gz | bcftools view -m2 -M2 -c1 -Oz -o snp_bcftools_FLT.vcf.gz -
bcftools stats snp_bcftools_FLT.vcf.gz > snp_bcftools_FLT.stats
tabix -f -p vcf snp_bcftools_FLT.vcf.gz
```

Variant annotation with [SnpEff v5.0](https://pcingola.github.io/SnpEff/se_buildingdb/)

```
snpEff_DIR=/home/software/snpEff
TMP=../DATA_temp
CONFIG=/home/software/snpEff/snpEff.config
# SnpEff
snpEff -c $CONFIG SC5314 snp_bcftools_FLT.vcf.gz > snp_bcftools_FLT_snpEff.vcf
bgzip -f snp_bcftools_FLT_snpEff.vcf
tabix -f -p vcf snp_bcftools_FLT_snpEff.vcf.gz
```

Splitting SNPs and INDELs
```
bcftools view -i '%TYPE="snp"' -Oz -o snp_bcftools_FLT_snpEff_SNP.vcf.gz snp_bcftools_FLT_snpEff.vcf.gz
bcftools view -i '%TYPE="indel"' -Oz -o snp_bcftools_FLT_snpEff_INDEL.vcf.gz snp_bcftools_FLT_snpEff.vcf.gz
```

Checking missingness of the data with [vcftools ](https://vcftools.github.io/man_latest.html)
```
VCF=snp_bcftools_FLT_snpEff.vcf.gz
vcftools --gzvcf $VCF --missing-indv --out snp_bcftools_FLT_snpEff
vcftools --gzvcf $VCF --missing-site --out snp_bcftools_FLT_snpEff
```


Additional example filters with vcftools:  
- min DP per sample = 10  
- min genotype quality = 20  
- max missing genotypes = 10%  
```
vcftools --gzvcf snp_bcftools_FLT_snpEff.vcf.gz --remove-filtered-all --minDP 10 --minGQ 20 --max-missing 0.9 --recode --recode-INFO-all --out snp_bcftools_FLT_snpEff_FLT2
```

### 9 SNP calling with GATK

If bams generate error - readjust RG field following [this site](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)  


Indexing reference
```
gatk CreateSequenceDictionary -R ref.fasta
```
Haplotype Caller (gatk-4.1.0.0)
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

### 10 SNP calling with Freebayes

Readjustment of RG fields in bam files using picard, if needed
```
picard AddOrReplaceReadGroups I=bams/${sample}_rmdup.bam O=bams/${sample}_RG.bam RGID=${sampleID} RGLB=${sampleID} RGPL=ILLUMINA RGPU=${flowcell.lane}.${sampleID} RGSM=${sampleID}
```
SNP calling with [freebayes](https://github.com/freebayes/freebayes) (v1.3.1-17-gaa2ace8)
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
tabix -p vcf freebayes_calls.f.vcf.gz
bcftools view -m2 -M2 --include 'TYPE="snp"' -O v -o freebayes_calls.f.snp.vcf freebayes_calls.f.vcf.gz
```

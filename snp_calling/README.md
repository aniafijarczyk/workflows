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


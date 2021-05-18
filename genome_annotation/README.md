### Gene annotation with Augustus


[1. Modelling repeats (RepeatModeler)](#1-Modelling-repeats)  
[2. Masking repeats (RepeatMasker)](#2-Masking-repeats)  
[3. *Ab initio* gene annotation (Augustus)](#3-Ab-initio-gene-annotation)  
[4. Generating RNAseq hints](#4-RNAseq-hints)
[5. Generating protein hints](#5-Protein-hints)
[6. Gene annotation with hints](#6-Gene-annotation-with-hints)

#### 1 Modelling repeats
Using RepeatModeler v2.0.1

```
genome=spades07_BM03_scaffolds.fasta
genome_name=BM03
mkdir -p RepMod_${genome_name}
cd RepMod_${genome_name}
BuildDatabase -name ${genome_name} ../${genome}
RepeatModeler -database ${genome_name} -pa 4 -LTRStruct
cd ..
```

#### 2 Masking repeats
Using RepeatMasker v4.1.0

```
genome=spades07_BM03_scaffolds.fasta
genome_name=BM03
/prg/RepeatMasker/4.1.0/util/queryRepeatDatabase.pl -species "fungi" > repeatmasker.fungi.fa
cd RepMod_${genome_name}
cat ${genome_name}-families.fa ../repeatmasker.fungi.fa > ${genome_name}.fulllib.fa
RepeatMasker -pa 4 -lib ${genome_name}.fulllib.fa ../${genome}
RepeatMasker -xsmall -pa 4 -lib ${genome_name}.fulllib.fa ../${genome}
cd ..
```

#### 3 *Ab initio* gene annotation
Using Augustus v3.3.2

```
genome=spades07_BM03_scaffolds_softmasked.fasta
genome_name=BM03
augustus --species=saccharomyces_cerevisiae_S288C ${genome} --outfile=augustus_abinitio_${genome_name}.gff --errfile=augustus_abinitio_${genome_name}.err --genemodel=partial --singlestrand=false
getAnnoFasta.pl augustus_abinitio_${genome_name}.gff --seqfile=${genome}
```
#### 4 RNAseq hints
```
strainID=BM03
genome=spades07_BM03_scaffolds_softmasked.fasta

# indexing [STAR v2.7.2b]
STAR --runMode genomeGenerate --genomeSAindexNbases 11 --genomeDir star_genome_${strainID} --genomeFastaFiles ${genome}

# mapping [STAR 2.7.2b & samtools v1.8]
STAR --genomeDir star_genome_${strainID} --readFilesIn ${fastq1},${fastq2},${fastq3} --readFilesCommand zcat --runThreadN 8
samtools view -Sb Aligned.out.sam > Aligned.out.bam
rm Aligned.out.sam
mv SJ.out.tab SJ.${strainID}.out.tab

# filtering primary alignment [bamtools v2.5.1]
bamtools filter -in Aligned.out.bam -out Aligned.out.f.bam -isPrimaryAlignment true

# sorting [samtools v1.8]
samtools sort -o Aligned.${strainID}.out.f.sort.bam Aligned.out.f.bam
samtools index Aligned.${strainID}.out.f.sort.bam
rm Aligned.out.f.bam

# generating exon hints [augustus v3.3.2, braker v2.1.2]
bam2wig Aligned.${strainID}.out.f.sort.bam  > rnaseq.${strainID}.wig
cat rnaseq.${strainID}.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 \
--minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.${strainID}.track --radius=4.5 --pri=4 --strand="." > rnaseq.${strainID}.gff

# generating intron hints [augustus v3.3.2, braker v2.1.2]
bam2hints --intronsonly --in=Aligned.${strainID}.out.f.sort.bam --out=introns.${strainID}.gff
filterIntronsFindStrand.pl $genome introns.${strainID}.gff --score > introns.${strainID}.f.gff
```

#### 5 Protein hints
```
genome=spades07_BM03_scaffolds_softmasked.fasta
strainID=BM03
proteins=yeast_proteins.fa
# alignment [braker v2.1.2]
startAlign.pl --genome=${genome} --prot=${proteins} --prg=gth --CPU=4

# getting prot hints [braker v2.1.2]
align2hints.pl --in=align_gth/gth.concat.aln --out=prot.hints.${strainID}.gff --prg=gth

# changing priority
less prot.hints.${strainID}.gff | sed 's/pri=4/pri=3/g' > prot.hints.pri3.${strainID}.gff
```

#### 6 Gene annotation with hints
```
cp */extrinsic/extrinsic.M.RM.E.W.P.cfg extrinsic.cfg
genome=spades07_BM03_scaffolds_softmasked.fasta
strainID=BM03
cat introns.${strainID}.f.gff prot.hints.pri3.${strainID}.gff > hints.${strainID}.gff

augustus --species=saccharomyces_cerevisiae_S288C \
  --extrinsicCfgFile=extrinsic.cfg --hintsfile=hints.${strainID}.gff \
  --allow_hinted_splicesites=atac \
  --softmasking=on ${genome} --outfile=augustus_hints_${strainID}.gff \
  --errfile=augustus_hints_${strainID}.err
```

#### 7 Joining gtf
```
species=BM03
# generating gtf files
cat augustus_abinitio_${species}.gff | grep -v "^#" > temp1.gtf
cat augustus_hints_${species}.gff | grep -v "^#" > temp2.gtf

# joining gene models
joingenes --genesets=temp1.gtf,temp2.gtf --priorities=1,2 --output=augustus_joined_${species}.gff --errordistance=878 -i
# "-i : If this flag is set, the program joines the stop_codons to the CDS"

# finding partial gene models and translate
getAnnoFasta.pl augustus_joined_${species}.gff --seqfile=${genome}
python findIncomplete.py
python translateCodingseq.py
```


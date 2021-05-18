### Gene annotation with Augustus


[1. Modelling repeats (RepeatModeler)](#1-Modelling-repeats)  
[2. Masking repeats (RepeatMasker)](#2-Masking-repeats)  
[3. *Ab initio* gene annotation (Augustus)](#3-Ab-initio-gene-annotation)  


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



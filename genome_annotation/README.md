### Gene annotation with Augustus


[1. Modelling repeats (RepeatModeler)](#1-Modelling-repeats)
[2. Masking repeats (RepeatMasker)](#2-Masking-repeats)
[3. *Ab initio* gene annotation (Augustus)](#3-Ab-initio-gene-annotation)


#### Modelling repeats
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

#### Masking repeats
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



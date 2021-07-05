### Scripts for manipulating vcf files



#### A scikit-allel wrapper to calculate dxy in windows and per chromosome/contig. To get an unbiased dxy, generate vcf file including non-variant sites.
```
python calculateDxy.py -i sample_files/sample.vcf.gz -p sample_files/pop1_mosaicbeer.tab -r sample_files/pop2_4eflv2.tab -w 25000 -s 25000
```


#### A scikit-allel wrapper to calculate distance of a sample/samples from the PCA centroid relative to other group of samples in non-overlapping windows of predetermined number of SNPs. Only sites with segreating biallelic variants and no missing genotypes are taken into account.
```
python calculatePCADist.py -i sample_files/sample.vcf.gz -p sample_files/pop1_mosaicbeer.tab -r sample_files/pop2_4eflv2.tab -w 100
```


### Scripts for manipulating vcf files



#### A scikit-allel wrapper to calculate dxy in windows and per chromosome/contig. To get an unbiased dxy, generate vcf file including non-variant sites.
```
python calculateDxy.py -i sample_files/sample.vcf.gz -p sample_files/pop1_mosaicbeer.tab -r sample_files/pop2_4eflv2.tab -w 25000 -s 25000
```

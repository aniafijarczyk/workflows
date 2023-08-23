#!/bin/bash

fastq=$1
nameis=$(basename $fastq)
zcat ${fastq} | head -n1 | awk -F":" -v name="${nameis}" '{print name,$3,$4,$10}'


#!/bin/bash

sample_id="Sample1"

fastqs=("../data/Sample1_L1_R1.fastq.gz" "../data/Sample1_L1_R2.fastq.gz" "../data/Sample1_L2_R1.fastq.gz" "../data/Sample1_L2_R2.fastq.gz") 

out="../data/${sample_id}_kb_out"

index="../ref/transcriptome.idx"
t2g="../ref/t2g.txt"

r1_len="$(zcat ${fastqs[0]} | head -80 | sed -n 2~4p | awk '{ print length }' | awk -F : '{sum+=$1} END {print sum/NR}' | cut -f1 -d".")"
if [[ $r1_len =~ 28 ]]
then
    tech="10xv3"
elif [[ $r1_len =~ 26 ]]
then
    tech="10xv2"
else
echo 'Unable to determine version'
fi

kb count -i $index -g $t2g -x $tech -o $out ${fastqs[@]}

Rscript ./reformat_raw.R ${sample_id}

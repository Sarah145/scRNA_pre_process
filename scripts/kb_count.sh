#!/bin/bash

### Replace "Sample1" below with your sample information in the two lines of code below  ###

# Assign a sample_id
sample_id="Sample1"

# Point to fastq files in R1/R2 pairs
fastqs=("../data/Sample_R1.fastq.gz" "../data/Sample1_L1_R2.fastq.gz" "../data/Sample1_L2_R1.fastq.gz" "../data/Sample1_L2_R2.fastq.gz") 

######

out="../data/${sample_id}_kb_out" # specify output directory

index="../ref/transcriptome.idx" # specify reference information
t2g="../ref/t2g.txt"

# detect chemistry version by checking the length of the first 10 R1s and taking the floor of the average
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

kb count -i $index -g $t2g -x $tech -o $out ${fastqs[@]} # run kallisto|bustools

Rscript ./reformat.R ${sample_id} # submit R script to reformat output

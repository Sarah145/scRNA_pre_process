#!/bin/bash

ref_dir="../ref/"
dna_fa="${ref_dir}Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf="${ref_dir}cellranger_genes.gtf.gz"

kb ref -i ${ref_dir}transcriptome.idx -g ${ref_dir}t2g.txt -f1 ${ref_dir}cdna.fa $dna_fa $gtf



#! /usr/bin/env Rscript

sample_id <- commandArgs(trailingOnly = TRUE)
print(paste('Working on', sample_id))

library(Matrix, quietly=T) # load libraries
library(DropletUtils, quietly=T)

raw_mtx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/matrix.mtx')) # load raw mtx
genes <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/genes.tsv'), sep = '\t', header = F) # load genes
rownames(raw_mtx) <- genes[,1] # attach gene_ids
colnames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/barcodes.tsv'), sep = '\t', header = F)[,1] # attach barcodes

out <- emptyDrops(raw_mtx) # get probability that each barcode is a cell
keep <- out$FDR <= 0.05 # define threshold probability for calling a cell
keep[is.na(keep)] <- FALSE
filt_mtx <- raw_mtx[,keep] # subset raw mtx to remove empty drops

write10xCounts(paste0('../data/', sample_id, '_kb_out/counts_filtered'), gene.symbol = genes[,2], filt_mtx, overwrite=T) # write out filtered results

print('Done!')

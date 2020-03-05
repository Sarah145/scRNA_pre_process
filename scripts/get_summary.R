#! /usr/bin/env Rscript

sample_id <- commandArgs(trailingOnly = TRUE)
print(paste('Working on', sample_id))

# load libraries
library(Matrix, quietly=T)
library(DropletUtils, quietly=T)
library(ggplot2, quietly=T)
library(scales, quietly=T)
library(rjson, quietly=T)
library(R2HTML, quietly=T)
source('./functions.R') # load bc_rank_plot and print_HTML functions

raw_mtx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/matrix.mtx')) # load raw mtx
filt_mtx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_filtered/matrix.mtx')) # load filtered mtx
kb_stats <- c(fromJSON(file = paste0('../data/', sample_id, '_kb_out/inspect.json')), 
		fromJSON(file = paste0('../data/', sample_id, '_kb_out/run_info.json'))) # load run info 
tech <- grep('10x(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) # determine chemistry version
seq_stats <- data.frame(stat = c('Sequencing technology', 'Number of reads processed', '% reads pseudoaligned', # get sequencing/alignment stats 
                                  '% reads valid after barcode correction'), 
                        value = prettyNum(c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
				  round(kb_stats$percentageReadsOnWhitelist,2)), big.mark = ','))
p_cnts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2) # calculate cell stats and save to df
med_cnts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated number of cells', '% counts in cells', 
                                 'Median counts per cell', 'Median genes per cell', 'Total genes detected'), 
                        value = prettyNum(c(ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
                                        med_genes_cell, tot_genes_detected), big.mark = ','))
stats <- barcodeRanks(raw_mtx) # get rank stats
raw_cells <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/barcodes.tsv'), header = F, sep ='\t')[,1] # load raw cells
filt_cells <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_filtered/barcodes.tsv'), header = F, sep ='\t')[,1] # load filtered cells
bc_rank_plot(stats = stats, raw_cells = raw_cells, filt_cells = filt_cells, save = paste0('../data/', sample_id, '_kb_out/barcode_rank.png')) # create barcode rank plot png
print_HTML(seq_stats = seq_stats, cell_stats = cell_stats, dir = paste0('../data/', sample_id, '_kb_out'), sample_id = sample_id) # output a HTML summary of the run

print('Done!')


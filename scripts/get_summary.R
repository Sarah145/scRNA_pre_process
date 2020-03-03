sample_id <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(DropletUtils)
library(ggplot2)
library(scales)
library(rjson)
library(R2HTML)
source('./functions.R')

raw_mtx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/matrix.mtx'))
filt_mx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_filtered/matrix.mtx'))
kb_stats <- c(fromJSON(file = paste0('../data/', sample_id, '_kb_out/inspect.json')), 
		fromJSON(file = paste0('../data/', sample_id, '_kb_out/run_info.json')))
tech <- grep('10x(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T)
seq_stats <- data.frame(stat = c('Sequencing technology', 'Number of reads processed', '% reads pseudoaligned', 
                                  '% reads valid after barcode correction'), 
                        value = prettyNum(c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
				  round(kb_stats$percentageReadsOnWhitelist,2)), big.mark = ','))
p_cnts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2)
med_cnts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated number of cells', '% counts in cells', 
                                 'Median counts per cell', 'Median genes per cell', 'Total genes detected'), 
                        value = prettyNum(c(ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
                                        med_genes_cell, tot_genes_detected), big.mark = ','))
stats <- barcodeRanks(raw_mtx)
raw_cells <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/barcodes.tsv'), header = F, sep ='\t')[,1]
filt_cells <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_filtered/barcodes.tsv'), header = F, sep ='\t')[,1]
bc_rank_plot(stats = stats, raw_cells = raw_cells, filt_cells = filt_cells, save = paste0('../data/', sample_id, '_kb_out/barcode_rank.png'))
printHTML(seq_stats = seq_stats, cell_stats = cell_stats, dir = paste0('../data/', sample_id, '_kb_out'), sample_id = sample_id)



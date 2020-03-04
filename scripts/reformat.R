sample_id <- commandArgs(trailingOnly = TRUE)

library(Matrix) # load libraries
library(DropletUtils)

raw_mtx <- as(t(readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.mtx'))), 'CsparseMatrix') # load mtx and transpose it
rownames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.genes.txt'), sep = '\t', header = F)[,1] # attach genes
colnames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.barcodes.txt'), header = F, sep = '\t')[,1] # attach barcodes

t2g <-  unique(read.csv('../ref/t2g.txt', sep = '\t', header=F)[,2:3]) # load t2g file
t2g <- data.frame(t2g[,2], row.names = t2g[,1])
gene_sym <- t2g[as.character(rownames(raw_mtx)),1] # get symbols for gene ids

write10xCounts(paste0('../data/', sample_id, '_kb_out/counts_unfiltered'), gene.symbol = gene_sym, raw_mtx, overwrite = T) # write results


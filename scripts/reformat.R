sample_id <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(DropletUtils)

raw_mtx <- as(t(readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.mtx'))), 'CsparseMatrix')
rownames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.genes.txt')), header = F, sep = '\t')[,1]
colnames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.barcodes.txt')), header = F, sep = '\t')[,1]
write10xCounts(paste0('../data/', sample_id, '_kb_out/counts_unfiltered'), raw_mtx, overwrite = T)


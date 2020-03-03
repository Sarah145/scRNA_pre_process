sample_id <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(DropletUtils)

raw_mtx <- as(t(readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.mtx'))), 'CsparseMatrix')
genes <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/genes.tsv'), sep = '\t', header = F)
t2g <-  unique(read.csv('../ref/t2g.txt', sep = '\t', header=F)[,2:3])
t2g <- data.frame(t2g[,2], row.names = t2g[,1])
rownames(raw_mtx) <- t2g[as.character(genes$V1),1]
colnames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/cells_x_genes.barcodes.txt')), header = F, sep = '\t')[,1]
write10xCounts(paste0('../data/', sample_id, '_kb_out/counts_unfiltered'), raw_mtx, overwrite = T)


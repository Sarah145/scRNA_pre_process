sample_id <- commandArgs(trailingOnly = TRUE)

library(Matrix)
library(DropletUtils)

raw_mtx <- readMM(paste0('../data/', sample_id, '_kb_out/counts_unfiltered/matrix.mtx'))
rownames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/genes.tsv'), sep = '\t', header = F)[,1]
colnames(raw_mtx) <- read.csv(paste0('../data/', sample_id, '_kb_out/barcodes.tsv'), sep = '\t', header = F)[,1]
out <- emptyDrops(raw_mtx)
keep <- out$FDR <= 0.05
keep[is.na(keep)] <- FALSE
filt_mtx <- raw_mtx[,keep]
write10xCounts(paste0('../data/', sample_id, '_kb_out/counts_filtered'), filt_mtx, overwrite=T) 

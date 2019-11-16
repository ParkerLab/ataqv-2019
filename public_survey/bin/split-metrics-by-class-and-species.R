#!/usr/bin/env Rscript

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
SPECIES <- args[2]
CLASS <- args[3]
METRICS_FILE <- args[4]
OUT <- args[5]

tmp <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
tmp <- tmp[tmp$organism == SPECIES,]
stopifnot(CLASS %in% c('single-cell', 'bulk'))
if (CLASS == 'single-cell') {
	tmp <- tmp[tmp$is_single_cell,]
} else {
	tmp <- tmp[!tmp$is_single_cell,]
}

metrics <- read.table(METRICS_FILE, head = F, as.is = T, sep = '\t')
metrics <- metrics[metrics$V1 %in% tmp$experiment,]

write.table(metrics, file=OUT, col.names = F, row.names = F, sep = '\t', quote = F, append = F)

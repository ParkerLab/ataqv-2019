#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
TSS_ENRICHMENT <- args[2]

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')
tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)
sample_info$library <- as.character(sample_info$library)

TSS_ENRICHMENT_COMPARISON <- read.table(TSS_ENRICHMENT, head = F, as.is = T, sep = '\t', col.names = c('file', 'metric', 'value')) %>%
  dplyr::mutate(library=gsub('(\\d+)-(.*).json.gz', '\\1', basename(file)),
                tss_list=gsub('(\\d+)-(.*).json.gz', '\\2', basename(file))) %>%
  dplyr::select(-file, -metric) %>%
  dplyr::rename(tss_enrichment=value)

limmin <- min(TSS_ENRICHMENT_COMPARISON$tss_enrichment)
limmax <- max(TSS_ENRICHMENT_COMPARISON$tss_enrichment)

p <- TSS_ENRICHMENT_COMPARISON %>%
  tidyr::spread(key = tss_list, value = tss_enrichment) %>%
  left_join(sample_info) %>%
  ggplot() + 
  geom_point(aes(x = all_refseq, y = housekeeping_refseq, color = tn5)) +
  theme_bw() +
  xlab('TSS enrichment using\nall RefSeq TSS') +
  ylab('TSS enrichment using RefSeq\nTSS of housekeeping genes') +
  scale_color_manual(values = tn5_colors, breaks = c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'), labels = c('0.2', '0.5', '0.66', '1', '1.5', '2', '5')) +
  guides(color = guide_legend(title = 'Relative [Tn5]')) +
  xlim(c(limmin, limmax)) +
  ylim(c(limmin, limmax)) +
  geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed')
pdf('tss_enrichment_housekeeping_vs_all.pdf', width = 6, height = 5)
print(p)
dev.off()

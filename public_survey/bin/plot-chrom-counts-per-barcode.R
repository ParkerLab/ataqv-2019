#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

args <- commandArgs(T)
COUNTS <- args[1]
SAMPLE_INFO <- args[2]

counts <- read.table(COUNTS, head = F, as.is = T, sep = '\t', col.names = c('barcode', 'chrom', 'count'), colClasses = c('character', 'character', 'numeric'))
counts <- counts %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(total_count_for_barcode=sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(proportion_for_chrom=count/total_count_for_barcode)

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t') %>%
  dplyr::filter(is_single_cell) %>%
  dplyr::select(experiment, organism, title) %>%
  dplyr::rename(barcode=experiment) %>%
  unique()
counts <- left_join(counts, sample_info)

expected_proportion_per_chromosome <- counts %>%
  dplyr::select(chrom, count, organism) %>%
  dplyr::group_by(chrom, organism) %>%
  dplyr::summarize(count=sum(count)) %>%
  dplyr::ungroup()
expected_proportion_per_chromosome$expectation <- expected_proportion_per_chromosome$count / sum(expected_proportion_per_chromosome$count)

# plot the maximum proportion of reads from any single chromosome...
sorted_chromosomes <- unique(paste0('chr', sort(as.numeric(gsub('^chr', '', counts$chrom)))))
counts$chrom <- factor(counts$chrom, levels=sorted_chromosomes, ordered = T)
max_proportion <- counts %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(max_proportion_from_a_single_chromosome=max(proportion_for_chrom)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(proportion_for_chrom==max_proportion_from_a_single_chromosome)
p <- ggplot(max_proportion) +
  geom_point(aes(x = total_count_for_barcode, y = max_proportion_from_a_single_chromosome, color=chrom), alpha = 0.2) +
  scale_x_log10() +
  scale_color_viridis_d() +
  ylab('Max fraction of reads from a single chromosome') +
  xlab('Read counts for barcode') +
  theme_bw() +
  facet_wrap(~organism)
pdf('max-proportion-from-a-single-chromosome.pdf', height = 5, width = 10)
p
dev.off()

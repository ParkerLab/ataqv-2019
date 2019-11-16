#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(ggrepel)
library(RColorBrewer)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
FRAGMENT_LENGTHS <- args[2]
TSS_COVERAGE <- args[3]

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
fragment_length_distributions <- read.table(FRAGMENT_LENGTHS, head = F, as.is = T, sep = '\t', col.names = c('experiment', 'fragment_length', 'count'))
fragment_length_distributions <- fragment_length_distributions %>%
	dplyr::group_by(experiment) %>%
	dplyr::mutate(fraction=count/sum(count)) %>%
	dplyr::ungroup()
tss_coverage <- read.table(TSS_COVERAGE, head = F, as.is = T, sep = '\t', col.names = c('experiment', 'position', 'coverage'))

fragment_length_distributions$experiment <- factor(fragment_length_distributions$experiment, levels = sample(unique(fragment_length_distributions$experiment)), ordered = T)
fragment_length_distributions <- fragment_length_distributions[order(fragment_length_distributions$experiment, fragment_length_distributions$fragment_length),]

tss_enrichment <- tss_coverage
tss_enrichment$experiment <- factor(tss_enrichment$experiment, levels = sample(unique(tss_enrichment$experiment)), ordered = T)
tss_enrichment <- tss_enrichment[order(tss_enrichment$experiment, tss_enrichment$position),]

# heterogeneity within studies
STUDIES_TO_SHOW <- c('PRJNA259243', 'PRJNA306754') # Qu et al, Ackermann
EXPERIMENTS_TO_SHOW <- sample_info$experiment[sample_info$project %in% STUDIES_TO_SHOW]
p1 <- fragment_length_distributions %>% dplyr::filter(experiment %in% EXPERIMENTS_TO_SHOW) %>% 
  left_join(sample_info %>% dplyr::select(experiment, project) %>% unique()) %>%
  ggplot() + geom_line(aes(x = fragment_length, y = fraction, color = experiment), alpha = 0.5) +
  guides(color=F) + xlab('Fragment length (bp)') + ylab('Fraction of reads') +
  coord_cartesian(xlim = c(0, 500)) +
  theme_bw() +
  facet_wrap(~project, ncol=1)

p2 <- tss_enrichment %>% dplyr::filter(experiment %in% EXPERIMENTS_TO_SHOW) %>% 
  left_join(sample_info %>% dplyr::select(experiment, project) %>% unique()) %>%
  ggplot() + geom_line(aes(x = position, y = coverage, color = experiment), alpha = 0.5) +
  guides(color=F) +
  xlab('Position relative to TSS (bp)') +
  ylab('Normalized read coverage') +
  theme_bw() +
  facet_wrap(~project, ncol=1)


pdf('public_intrastudy_heterogeneity_fragment_length_distributions.pdf', height = 4, width = 3)
print(p1)
dev.off()

pdf('public_intrastudy_heterogeneity_tss_enrichment.pdf', height = 4, width = 3)
print(p2)
dev.off()


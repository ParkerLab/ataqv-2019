#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(latex2exp)

args <- commandArgs(T)
PC_SCORES <- args[1]
METRICS <- args[2] 
PREFIX <- args[3]
PROJECT <- args[4]
PROJECT_COLORS <- args[5]

metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric', 'value')) %>%
  dplyr::mutate(value=as.numeric(value))
pcs <- read.table(PC_SCORES, head = T, as.is = T, sep = '\t')
metrics <- left_join(pcs, metrics)

project_colors <- read.table(PROJECT_COLORS, head = T, as.is = T, sep = '\t', comment.char = '')
COLOR <- project_colors$color[project_colors$project==PROJECT]

# remove metrics that are invariant
invariant <- metrics %>%
  dplyr::group_by(metric) %>%
  dplyr::summarize(unique_vals=length(unique(value))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(unique_vals==1)
metrics <- metrics[!metrics$metric %in% invariant$metric,]
correlations <- metrics %>%
  dplyr::group_by(metric) %>%
  dplyr::summarize(abs_pearsons_cor=cor(value, PC1)^2)
correlations <- correlations[order(correlations$abs_pearsons_cor),]
correlations$metric <- factor(correlations$metric, levels = correlations$metric, ordered = T)

p <- ggplot(correlations) +
  geom_bar(aes(x = metric, y = abs_pearsons_cor), stat = 'identity') +
  theme_bw() +
  coord_flip() +
  xlab('') +
  ylab(TeX("$r^2$ between metric and PC1")) +
  ggtitle(glue('In project ID {PROJECT}')) +
  theme(plot.title = element_text(hjust = 0.5))
pdf(glue('{PREFIX}cor_barplot.pdf'), height = 6, width = 5)
p
dev.off()


p <- ggplot(metrics) +
  geom_point(aes(x = PC1, y = value), alpha = 0.4, color = COLOR) +
  facet_wrap(~metric, scales = 'free') +
  theme_bw() +
  xlab('PC1 score') +
  ylab('QC metric value') +
  ggtitle(glue('In project ID {PROJECT}')) +
  theme(plot.title = element_text(hjust = 0.5))
pdf(glue('{PREFIX}all_scatter.pdf'), height = 15, width = 15)
p
dev.off()

tss <- metrics[metrics$metric=='tss_enrichment',]

p <- ggplot(tss) +
  geom_point(aes(x = PC1, y = value), alpha = 0.4, color = COLOR) +
  theme_bw() +
  xlab('PC1 score') +
  ylab('TSS enrichment') +
  ggtitle(glue('In project ID {PROJECT}')) +
  theme(plot.title = element_text(hjust = 0.5))
pdf(glue('{PREFIX}tss_enrichment_scatter.pdf'), height = 2, width = 2)
p
dev.off()



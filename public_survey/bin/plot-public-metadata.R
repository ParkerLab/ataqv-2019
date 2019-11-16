#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

# need sample info, metrics, fragment lengths, tss coverage
args <- commandArgs(T)
SAMPLE_INFO <- args[1]
PROJECT_COLORS <- args[2]

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
project_colors <- read.table(PROJECT_COLORS, head = T, as.is = T, sep = '\t', comment.char = '')
COLORS <- project_colors$color
names(COLORS) <- project_colors$project

# characterize public metadata for Fig 1
tmp <- sample_info %>% dplyr::group_by(experiment, project, organism, is_single_cell) %>%
  dplyr::summarize(n_reads=sum(as.numeric(spots)))
tmp$is_single_cell <- ifelse(tmp$is_single_cell, 'single cell', 'bulk')
tmp$Species <- ifelse(tmp$organism=='Mus musculus','mouse','human')
tmp$label <- paste(tmp$Species, tmp$is_single_cell)
counts_per_group <- as.data.frame(table(tmp$label))
colnames(counts_per_group) <- c('label', 'count')
tmp <- left_join(tmp, counts_per_group)
tmp$label <- paste0(tmp$label, '\n(N = ', tmp$count, ')')
label_levels_df <- unique(tmp[,c('label', 'is_single_cell', 'Species')])
label_levels_df <- label_levels_df[rev(order(label_levels_df$is_single_cell, label_levels_df$Species)),]
label_levels_df <- label_levels_df[c(2, 1, 4, 3),]
tmp$label <- factor(tmp$label, levels = label_levels_df$label, ordered = T)
p <- ggplot(tmp) + geom_jitter(aes(x = label, y = log10(n_reads), color = project), stat = 'identity', width = 0.4, height = 0, alpha = 0.3, stroke = 0) +
  scale_color_manual(values = COLORS, guide = F) +
  xlab('') +
  ylab('Log10(read pairs in library)') +
  theme_bw() +
  coord_flip()
pdf('public_meta.pdf', height = 2.5, width = 3.5)
p
dev.off()

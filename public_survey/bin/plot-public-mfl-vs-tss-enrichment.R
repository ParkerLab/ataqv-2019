#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
METRICS <- args[2]
PROJECT_COLORS <- args[3]

project_colors <- read.table(PROJECT_COLORS, head = T, as.is = T, sep = '\t', comment.char='')
COLORS <- project_colors$color
names(COLORS) <- project_colors$project

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('experiment', 'metric', 'value'))

tmp <- metrics %>%
  tidyr::spread(key = metric, value = value, convert = T) %>%
  left_join(sample_info %>% dplyr::select(experiment, project)) %>%
  unique()
p <- ggplot(tmp) + geom_point(aes(x = median_fragment_length, y = tss_enrichment, color = project), alpha = 0.5, stroke = 0) +
  geom_text(aes(x = 250, y = 20, label = paste0('N = ', nrow(tmp))), data = data.frame()) +
  theme_bw() + xlab('Median fragment length (bp)') + ylab('TSS enrichment') +
  scale_color_manual(values = COLORS, guide = F) +
  theme(axis.title.x = element_text(color = rgb(0, 174, 239, maxColorValue = 255)),
        axis.title.y = element_text(color = rgb(139, 94, 60, maxColorValue = 255)))
pdf('public_median_fragment_length_vs_tss_enrichment.pdf', height = 3, width = 3)
print(p)
dev.off()

total_number_public_read_pairs <- round(sum(sample_info$spots)/1e9, 1)
number_bulk_samples_shown <- length(unique(sample_info$experiment[!sample_info$is_single_cell]))

FOLD_DIFFERENCES <- metrics %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(value=as.numeric(value)) %>%
  dplyr::summarize(max_value=max(value),
                   min_value=min(value),
                   fc=max_value/min_value)
  
# TODO: update this in the manuscript
print("FOR MANUSCRIPT: FOLD DIFFERENCES IN PUBLIC DATA")
print(FOLD_DIFFERENCES)

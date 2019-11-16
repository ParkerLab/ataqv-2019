#!/usr/bin/env Rscript
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(ggrepel)
library(RColorBrewer)
library(broom)

tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
OVERLAP <- args[2]
TOTAL_READ_COUNTS <- args[3]

total_read_counts <- read.table(TOTAL_READ_COUNTS, head = F, as.is = T, sep = '\t', col.names = c('total_reads', 'library'), colClasses = c('numeric', 'character'))

# tn5
sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)
sample_info$library <- as.character(sample_info$library)
sample_info <- unique(sample_info[,c('library', 'tn5')])


overlap <- read.table(OVERLAP, head = F, as.is = T, sep = '\t', col.names = c('library', 'experiment', 'count'), colClasses = c('character', 'character', 'numeric')) %>%
	dplyr::mutate(factor=gsub('\\..*', '', experiment), 
		      accession=gsub('.*\\.', '', experiment)) %>%
	left_join(total_read_counts) %>%
	left_join(sample_info) %>%
	dplyr::mutate(proportion_of_reads_overlapping_peak=count/total_reads) %>%
	dplyr::group_by(factor, accession, tn5) %>%
	dplyr::summarize(prop=mean(proportion_of_reads_overlapping_peak))


prop_at_1x <- overlap %>% dplyr::filter(tn5=='1X') %>% dplyr::select(-tn5) %>% dplyr::rename(base_prop=prop)
tmp <- left_join(overlap, prop_at_1x) %>% mutate(normalized_prop=prop/base_prop)
tmp$label <- paste0(tmp$factor, ' (', tmp$accession, ')')
p <- ggplot(tmp %>% dplyr::filter(tn5!='1X')) + geom_tile(aes(x = gsub('X', '', tn5), y = factor, fill = normalized_prop)) +
  scale_fill_gradient2(low = 'blue', mid = 'white', 'high' = 'red', midpoint = 1) +
  xlab('Relative [Tn5]') + ylab('') +
  guides(fill = guide_colorbar(title = 'Proportion of reads overlapping ChIP-seq peaks\n(relative to relative [Tn5] = 1)')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = 'bottom') +
  coord_flip()
pdf('tn5_tf_overlap.pdf', height = 4, width = 8)
p
dev.off()

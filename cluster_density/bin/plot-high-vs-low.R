#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
METRICS <- args[2]

tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')

metric_names <- c(hqaa = 'High-quality autosomal reads',
                  percent_autosomal_duplicate = '% autosomal duplicate',
                  percent_mitochondrial = '% mitochondrial reads',
                  median_fragment_length = 'Median fragment length (bp)',
                  total_peaks = 'Number of peaks',
                  total_reads = 'Total number of reads',
                  hqaa_overlapping_peaks_percent = '% high-quality autosomal\nreads overlapping peaks',
                  tss_enrichment = 'TSS enrichment',
                  short_mononucleosomal_ratio = 'Short : mononucleosomal reads',
                  fragment_length_distance = 'Fragment length distance')


sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)
sample_info$library <- as.character(sample_info$library)
metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric', 'value'))

# global functions
# for converting between numeric Tn5 and Factor Tn5
# Example usage:
# as.factor.tn5('2')
# as.numeric.tn5(as.factor.tn5('2'))
tn5_conversions <- data.frame(num=c(0.2, 0.5, 0.66, 1, 1.5, 2, 5))
tn5_conversions$fac <- factor(paste(tn5_conversions$num, 'X', sep=''), paste(tn5_conversions$num, 'X', sep = ''), ordered = T)
rownames(tn5_conversions) <- as.character(tn5_conversions$fac)

as.numeric.tn5 <- function(tn5) {
  return(as.numeric(as.character(gsub('X', '', tn5))))
}

as.factor.tn5 <- function(tn5) {
  if (is.character(tn5)) {
    return(tn5_conversions[tn5,'fac'])
  } else if (is.numeric(tn5)) {
    return(tn5_conversions[paste(tn5, 'X', sep=''),'fac'])
  }
}

# cluster density plots
metrics.high_vs_low <- metrics %>%
  tidyr::separate(col = library, into = c('library', 'run')) %>%
  dplyr::mutate(run=ifelse(run==1830, 'low', 'high')) %>% 
  tidyr::spread(key = run, value = value)

metrics.high_vs_low %>%
  dplyr::mutate(difference=high-low) %>%
  dplyr::group_by(metric) %>%
  dplyr::summarize(avg_diff=mean(difference))


# set the color scheme
colors <- tn5_colors
names(colors) <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

for(each_metric in unique(metrics.high_vs_low$metric)) {
  tmp <- metrics.high_vs_low %>% dplyr::filter(metric==each_metric) %>%
    left_join(sample_info %>% dplyr::select(library, replicate, tn5) %>% dplyr::mutate(library=gsub('___.*', '', library)) %>% unique())
  
  xlabel <- paste(metric_names[each_metric], '\n(low cluster density)')
  ylabel <- paste(metric_names[each_metric], '\n(high cluster density)')
  lim.max <- max(c(tmp[,'low'], tmp[,'high']))
  lim.min <- min(c(tmp[,'low'], tmp[,'high']))
  p <- ggplot(tmp) + geom_point(aes(x = low, y = high, color = tn5, shape = replicate)) +
    theme_bw() +
    xlab(xlabel) +
    ylab(ylabel) +
    coord_cartesian(xlim=c(lim.min, lim.max),
                    ylim=c(lim.min, lim.max)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    scale_color_manual(values = colors)
  
  OUT <- glue('cd_high_vs_low_{each_metric}.pdf')
  pdf(OUT, height = 4, width = 5)
  print(p)
  dev.off()
}

# TODO:
# quantifications in the paper
# "The average difference in median fragment length between sequencing runs was XXX bps"
tmp <- metrics.high_vs_low %>% dplyr::filter(metric=='median_fragment_length') %>% dplyr::mutate(diff=low-high)
mean(tmp$diff)


# TSS enrichment was consistently higher in the high cluster density sequencing run (average difference of XXX
tmp <- metrics.high_vs_low %>% dplyr::filter(metric=='tss_enrichment') %>% dplyr::mutate(diff=high-low)
mean(tmp$diff)

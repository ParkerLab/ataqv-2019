#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("--sample_info"), action = 'store', type = 'character', help = '[Required] Path to sample info file'),
  make_option(c("--fld"), action = 'store', type = 'character', help = '[Required] Path to file of fragment length counts (library, length, count)'),
  make_option(c("--tss_coverage"), action = 'store', type = 'character', help = '[Required] Path to file of TSS coverage')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# GLOBALS -- should be shared across plotting scripts
tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')

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

# END GLOBALS

sample_info <- read.table(opts$sample_info, head = T, as.is = T, sep = '\t') %>%
	dplyr::mutate(library=as.character(library))
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)

fragment_lengths <- read.table(opts$fld, head = F, as.is = T, sep = '\t', col.names = c('library', 'fragment_length', 'count'), colClasses=c('character', 'numeric', 'numeric')) %>% dplyr::group_by(library) %>% dplyr::mutate(fraction=count/sum(count)) %>% dplyr::ungroup()
tss_coverages <- read.table(opts$tss_coverage, head = F, as.is = T, sep = '\t', col.names = c('library', 'position', 'coverage'), colClasses=c('character', 'numeric', 'numeric'))


fragment_lengths <- left_join(fragment_lengths, sample_info %>% dplyr::select(library, tn5, replicate))
tmp <- fragment_lengths %>% dplyr::filter(replicate=='rep1') %>%
  dplyr::filter(tn5 %in% c('0.5X', '1X', '2X'))
tmp$tn5 <- factor(tmp$tn5, levels = c('2X', '1X', '0.5X'), ordered = T)
fragment_length_plot <- ggplot(tmp) +
  geom_line(aes(x = fragment_length, y = fraction, color = tn5)) +
  scale_color_manual(values = c('2X' = tn5_colors[6], '0.5X' = tn5_colors[2], '1X' = tn5_colors[4]), labels = c('2', '1', '0.5')) +
  theme_bw() +
  coord_cartesian(xlim=c(0, 500)) +
  guides(color = guide_legend(title = 'Relative\n[Tn5]', override.aes = list(shape = 22))) +
  xlab("Fragment length (bp)") + ylab("Fraction of fragments") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=10), legend.title=element_text(size=12),
        legend.position = c(0.7, 0.7))

pdf('tn5_fragment_length_distributions.pdf', height = 3, width = 3)
print(fragment_length_plot)
dev.off()



tss_coverages <- left_join(tss_coverages, sample_info %>% dplyr::select(library, tn5, replicate))
tmp <- tss_coverages %>% dplyr::filter(replicate=='rep1') %>%
  dplyr::filter(tn5 %in% c('0.5X', '1X', '2X'))
tmp$tn5 <- factor(tmp$tn5, levels = c('2X', '1X', '0.5X'), ordered = T)
tss_coverages_plot <- ggplot(tmp) +
  geom_line(aes(x = position, y = coverage, color = tn5)) +
  scale_color_manual(values = c('0.5X' = tn5_colors[2], '1X' = tn5_colors[4], '2X' = tn5_colors[6]), labels = c('2', '1', '0.5')) +
  theme_bw() +
  xlab("Position relative to TSS") + ylab("Normalized read coverage") +
  ylim(c(0, max(tss_coverages$coverage))) +
  guides(color = guide_legend(title = 'Relative\n[Tn5]')) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.position = c(0.8, 0.7), legend.text = element_text(size=10), legend.title=element_text(size=12))

pdf('tn5_tss_enrichment.pdf', height = 3, width = 3)
print(tss_coverages_plot)
dev.off()


#!/usr/bin/env Rscript
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

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

# Compare TSS enrichment using different methods
# calculate the baseline for each method
args <- commandArgs(T)
SAMPLE_INFO <- args[1]
COVERAGE_FILES <- args[2:length(args)]

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$library <- as.character(sample_info$library)
sample_info$tn5 <- as.factor.tn5(sample_info$tn5)

load_coverage_file <- function(f) {
	lib <- gsub('(.*)\\.(.*)\\.txt', '\\1', basename(f))
	method <- gsub('(.*)\\.(.*)\\.txt', '\\2', basename(f))
	tmp <- read.table(f, head = F, as.is = T, sep = '\t', col.names = c('position', 'coverage'), colClasses = c('numeric', 'numeric'))
	tmp$library <- lib
	tmp$method <- method
	return(tmp)
}

tss_coverage <- bind_rows(lapply(COVERAGE_FILES, load_coverage_file))
tss_coverage$method[tss_coverage$method=='ataqv'] <- 'Coverage using fragment'
tss_coverage$method[tss_coverage$method=='cutsite'] <- 'Coverage using cutsite +/- 1 bp'
tss_coverage$method[tss_coverage$method=='encode'] <- 'Coverage using cutsite +/- 1/2 read length'

# First, just make a coverage plot for one library
SHOW_LIBRARY <- sample_info$library[sample_info$tn5=='1X' & sample_info$replicate=='rep1']
p <- tss_coverage %>%
  dplyr::filter(library==SHOW_LIBRARY) %>%
  ggplot() +
  geom_line(aes(x = position, y = coverage, color = method)) +
  scale_color_viridis_d() +
  theme_bw() +
  xlab('Position relative to TSS (bp)') +
  ylab('Normalized coverage') +
  guides(color = guide_legend(title=''))
pdf('tss_coverage_different_methods.pdf', width = 8, height = 3.5)
print(p)
dev.off()

# show the positional stability under each method
p <- tss_coverage %>%
  dplyr::group_by(library, method) %>%
  dplyr::mutate(max_coverage=max(coverage)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(coverage==max_coverage) %>%
  left_join(sample_info) %>%
  dplyr::mutate(label=glue('{tn5} Tn5, {replicate}')) %>%
  ggplot() +
  geom_point(aes(x = position, y=label, color = method)) +
  theme_bw() +
  scale_color_viridis_d() +
  xlab('Position of max coverage relative to TSS') +
  ylab('Library') +
  guides(color=F)
pdf('tss_coverage_positional_stability.pdf', width = 5, height = 8)
print(p)
dev.off()


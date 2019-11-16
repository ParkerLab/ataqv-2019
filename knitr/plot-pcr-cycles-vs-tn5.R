#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(ggrepel)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]

tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')
tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

# load in key data files/set key universal variables
# cluster density
sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)
sample_info$library <- as.character(sample_info$library)

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

# plot number of pcr cycles
tmp <- sample_info %>% dplyr::select(tn5, replicate, pcr_cycles) %>% unique()
tmp$tn5 <- log2(as.numeric.tn5(tmp$tn5))
p <- ggplot(tmp) + 
  geom_line(aes(x = tn5, y = pcr_cycles, color = gsub('rep', '', replicate))) +
  geom_point(aes(x = tn5, y = pcr_cycles, color = gsub('rep', '', replicate))) +
  scale_color_manual(values=tn5_replicate_colors) +
  theme_bw() + xlab('log2(Relative [ Tn5 ])') + ylab('Number of PCR cycles for library') +
  guides(color=guide_legend(title='Replicate'))
pdf('cd_pcr_cycles.pdf', width = 6, height = 5)
print(p)
dev.off()

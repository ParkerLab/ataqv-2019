#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--sample-info"), dest = 'sample_info', action = 'store', type = 'character', help = '[Required] Path to sample info file'),
  make_option(c("--rda"), action = 'store', type = 'character', help = '[Required] Path to pca.Rda'),
  make_option(c("--out"), action = 'store', type = 'character', help = '[Required] PDF to be written')
)

description <- 'Plot PCA'

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = description)
opts <- parse_args(option_parser)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


# GLOBALS
tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')
# END GLOBALS



sample_info <- read.table(opts$sample_info, head = T, as.is = T, sep = '\t')
sample_info$tn5 <- factor(sample_info$tn5, levels = tn5_levels, ordered = T)
sample_info$library <- as.character(sample_info$library)

load(opts$rda)
rpkm.prcomp <- pca
pca <- as.data.frame(rpkm.prcomp$x)

variances <- rpkm.prcomp$sdev^2 / sum(rpkm.prcomp$sdev^2)

sample_info.tmp <- sample_info %>% dplyr::select(-description) %>% tidyr::unite(col = info, library, tn5, experiment, replicate, sep = '::', remove = F)
sample_info.tmp$replicate <- factor(sample_info.tmp$replicate, ordered = F)
sample_info.tmp$day <- ifelse(sample_info.tmp$experiment == 'ex1', '1', '2')

pca$library <- rownames(pca)
pca <- left_join(pca, sample_info.tmp)

make_axis_label <- function(pc) {
  # pc should be an integer
  percentage <- round(variances[pc]*100, 2)
  percentage_string <- paste('(', percentage, '% of variance)', sep = '')
  return(paste('PC', pc, percentage_string))
}

p <- ggplot(pca) + geom_point(aes(x = PC1, y = PC2, color = tn5, shape = replicate, size = day)) + theme_bw() +
  xlab(make_axis_label(1)) +
  ylab(make_axis_label(2)) + theme(legend.box = 'horizontal') +
  scale_color_manual(values = tn5_colors, labels = c('0.2', '0.5', '0.66', '1', '1.5', '2', '5')) +
  scale_size_manual(values = c('1' = 1, '2' = 5)) +
  guides(color = guide_legend(title = 'Relative [Tn5]'), size = guide_legend(title = 'Day'), shape = guide_legend(title = 'Replicate'))

pdf('tn5_pca.pdf', width = 6, height = 3)
print(p)
dev.off()

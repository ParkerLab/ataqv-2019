#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("--sample_info"), action = 'store', type = 'character', help = '[Required] Path to sample info file'),
  make_option(c("--metrics"), action = 'store', type = 'character', help = '[Required] Path to file of ataqv metrics (library, metric, value)'),
  make_option(c("--plot_metric"), action = 'store', type = 'character', help = '[Required] Plot this metric'),
  make_option(c("--ylab"), action = 'store', type = 'character', help = '[Optional] Y axis label'),
  make_option(c("--out"), action = 'store', type = 'character', help = '[Required] Name of output file')
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

metrics <- read.table(opts$metrics, head = F, as.is = T, sep = '\t', col.names=c('library', 'metric', 'value'), colClasses=c('character', 'character', 'numeric')) %>%
	dplyr::filter(metric==opts$plot_metric) %>%
	left_join(sample_info)
metrics$tn5 <- log2(as.numeric.tn5(metrics$tn5))
metrics$replicate <- gsub('rep', '', metrics$replicate)

p <- ggplot(metrics) +
  geom_point(aes(x = tn5, y = value, color = replicate)) +
  geom_line(aes(x = tn5, y = value, color = replicate)) +
  scale_color_manual(values=tn5_replicate_colors) +
  xlab("log2(Relative [Tn5])") +
  ylab(opts$ylab) +
  theme_bw() +
  theme(axis.text=element_text(size=13))+
  theme(axis.title=element_text(size=12)) +
  guides(color = guide_legend(ncol = 1, title = element_text('Replicate')))

pdf(opts$out, height = 3, width = 4)
print(p)
dev.off()

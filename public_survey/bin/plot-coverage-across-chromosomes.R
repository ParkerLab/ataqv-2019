#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(T)
CHROM_COVERAGE <- args[1]
OUT <- args[2]

coverage <- read.table(CHROM_COVERAGE, head = F, as.is = T, sep = '\t', col.names = c('chrom', 'start', 'end', 'count'))

is_autosome <- function(chrom) {
  return(grepl('^chr\\d+$', chrom))
}

coverage <- coverage[is_autosome(coverage$chrom),]
sorted_chromosomes <- paste0('chr', sort(as.numeric(gsub('^chr', '', unique(coverage$chrom)))))
coverage$chrom <- factor(coverage$chrom, levels = sorted_chromosomes, ordered = T)
coverage$col <- ifelse(as.numeric(gsub('chr', '', coverage$chrom)) %% 2 == 0, 'even', 'odd')
coverage <- coverage[order(coverage$chrom, coverage$start),]
coverage$index <- 1:nrow(coverage)

labels <- coverage %>%
  dplyr::select(chrom, start, end) %>%
  dplyr::group_by(chrom) %>%
  dplyr::summarize(bins=n(),
                   middle=bins/2) %>%
  dplyr::ungroup()
labels$label_position <- c(0, cumsum(labels$bins[1:(nrow(labels)-1)])) + labels$middle


p <- ggplot(coverage) +
  geom_point(aes(x = index, y = count, color=col)) +
  # geom_line(aes(x = index, y = count)) +
  scale_x_continuous(breaks=labels$label_position, labels = labels$chrom) +
  scale_color_manual(values=c('odd' = 'blue', 'even' = 'red'), guide=F) +
  xlab('') +
  ylab('Reads in genomic bin') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))
pdf(OUT, height = 4, width = 10)
print(p)
dev.off()


#p <- ggplot(coverage) +
#  geom_density(aes(x = log(count+1), fill=chrom), alpha=0.1)
#p

#p <- ggplot(coverage) +
#  geom_boxplot(aes(x = chrom, y = log(count+1)))
#p


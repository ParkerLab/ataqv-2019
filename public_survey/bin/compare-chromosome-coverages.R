#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

args <- commandArgs(T)
NONOUTLIERS <- args[1]
CHROM_COVERAGE <- args[2:length(args)]
OUTLIER <- 'SRX859231'

nonoutliers <- read.table(NONOUTLIERS, head = F, as.is = T, sep = '\t')[,1]

read_coverage_file <- function(f) {
  tmp <- read.table(f, head = F, as.is = T, sep = '\t', col.names = c('chrom', 'start', 'end', 'count'))
  tmp$library <- gsub('.*(SRX\\d+).*', '\\1', basename(f))
  return(tmp)
}

coverage <- bind_rows(lapply(CHROM_COVERAGE, read_coverage_file))
coverage <- coverage[coverage$library %in% c(nonoutliers, OUTLIER),]

is_autosome <- function(chrom) {
  return(grepl('^chr\\d+$', chrom))
}

coverage <- coverage[is_autosome(coverage$chrom),]
sorted_chromosomes <- paste0('chr', sort(as.numeric(gsub('^chr', '', unique(coverage$chrom)))))
coverage$chrom <- factor(coverage$chrom, levels = sorted_chromosomes, ordered = T)
coverage$col <- ifelse(as.numeric(gsub('chr', '', coverage$chrom)) %% 2 == 0, 'even', 'odd')
coverage <- coverage[order(coverage$chrom, coverage$start),]
coverage <- coverage %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(index=1:length(count)) %>%
  dplyr::ungroup()

labels <- coverage %>%
  dplyr::select(chrom, start, end) %>%
  unique() %>%
  dplyr::group_by(chrom) %>%
  dplyr::summarize(bins=n(),
                   middle=bins/2) %>%
  dplyr::ungroup()
labels$label_position <- c(0, cumsum(labels$bins[1:(nrow(labels)-1)])) + labels$middle
labels$chrom <- gsub('chr', '', labels$chrom)

coverage$smoothed <- NA

# normalize to adjusted for differences in sequencing depth...
normalization_factor <- coverage %>%
  dplyr::filter(end%%coverage$end[1]==0) %>%
  dplyr::group_by(library) %>%
  dplyr::summarize(nf=sum(count))
coverage <- left_join(coverage, normalization_factor) %>%
  dplyr::mutate(normalized_coverage=count/nf)

coverage$smoothed <- NA
for(i in unique(coverage$chrom)) {
  for(l in unique(coverage$library)) {
    coverage$smoothed[coverage$chrom==i & coverage$library==l] <- loess(normalized_coverage ~ index, span=0.2, degree=1, normalized=F, data = coverage[coverage$chrom==i & coverage$library==l,])$fitted
  }
}


colors <- c(rep('grey', length(nonoutliers)), 'red')
names(colors) <- c(nonoutliers, OUTLIER)
alphas <- c(rep(0.1, length(nonoutliers)), 1)
names(alphas) <- c(nonoutliers, OUTLIER)

chromosome_boundaries <- unique(coverage[coverage$start==0,'index'])

p <- ggplot(coverage) +
  geom_line(aes(x = index, y = normalized_coverage, color=library, group=library, alpha=library), data=coverage[coverage$library==OUTLIER,]) +
  geom_line(aes(x = index, y = normalized_coverage, color=library, group=library, alpha=library), data=coverage[coverage$library!=OUTLIER,]) +
  scale_x_continuous(breaks=labels$label_position, labels = labels$chrom) +
  scale_color_manual(values=colors, guide=F) +
  scale_alpha_manual(values=alphas, guide=F) +
  xlab('Chrom.') +
  ylab('Norm. coverage\nof genomic bin') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  geom_vline(aes(xintercept = index), color='black', linetype='dashed', data=chromosome_boundaries) +
  coord_cartesian(ylim=c(0, 0.012))
pdf('compare-chromosome-coverages.pdf', height = 2, width = 5)
print(p)
dev.off()

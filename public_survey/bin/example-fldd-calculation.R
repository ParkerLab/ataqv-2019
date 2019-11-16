#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)

# example FLD calculation
args <- commandArgs(T)
FRAGMENT_LENGTHS <- args[1]

fragment_length_distributions <- read.table(FRAGMENT_LENGTHS, head = F, as.is = T, sep = '\t', col.names = c('experiment', 'fragment_length', 'count')) %>%
	dplyr::group_by(experiment) %>%
	dplyr::mutate(fraction=count/sum(count)) %>%
	dplyr::ungroup()

undertransposed <- 'SRX298004'
overtransposed <- 'SRX298009'
reference <- 'SRX298000'

undertransposed <- fragment_length_distributions %>% dplyr::filter(experiment==undertransposed) %>% dplyr::select(fragment_length, count, fraction)
undertransposed$cumulative <- cumsum(undertransposed$fraction)
undertransposed$sample <- 'undertransposed'

overtransposed <- fragment_length_distributions %>% dplyr::filter(experiment==overtransposed) %>% dplyr::select(fragment_length, count, fraction)
overtransposed$cumulative <- cumsum(overtransposed$fraction)
overtransposed$sample <- 'overtransposed'

reference <- fragment_length_distributions %>% dplyr::filter(experiment==reference) %>% dplyr::select(fragment_length, count, fraction)
reference$cumulative <- cumsum(reference$fraction)
reference$sample <- 'reference'

all <- rbind(overtransposed, undertransposed, reference)
all <- tidyr::gather(all, key = stat, value = value, fraction:cumulative)

# Create the vertical dashed lines representing the statistic
diff.o <- abs(reference$cumulative - overtransposed$cumulative)
point_of_greater_difference <- reference$fragment_length[diff.o==max(diff.o)]
ymax <- overtransposed$cumulative[overtransposed$fragment_length==point_of_greater_difference]
ymin <- reference$cumulative[reference$fragment_length==point_of_greater_difference]
line_between_reference_and_overtransposed <- data.frame(x = point_of_greater_difference, ymax = ymax, ymin = ymin, stat = 'cumulative')

diff.o <- abs(reference$cumulative - undertransposed$cumulative)
point_of_greater_difference <- reference$fragment_length[diff.o==max(diff.o)]
ymax <- undertransposed$cumulative[undertransposed$fragment_length==point_of_greater_difference]
ymin <- reference$cumulative[reference$fragment_length==point_of_greater_difference]
line_between_reference_and_undertransposed <- data.frame(x = point_of_greater_difference, ymax = ymax, ymin = ymin, stat = 'cumulative')

colors <- c('#66c2a5','#fc8d62','#8da0cb')
names(colors) <- rev(c('overtransposed', 'reference', 'undertransposed'))

# now piece together the final plot
cumulative <- ggplot(all[all$stat=="cumulative",]) + geom_step(aes(x = fragment_length, y = value, color = sample)) + theme_bw() + xlab("Fragment length (bp)") +
  coord_cartesian(xlim=c(0, 500)) +
  ylab("F(fragment length)") +
  guides(color = guide_legend(title = '')) +
  scale_color_manual(values = colors) +
  theme(legend.position = c(0.75, 0.35), legend.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
  geom_linerange(aes(x = x, ymax = ymax, ymin = ymin), data = line_between_reference_and_overtransposed, linetype = 'dashed') +
  geom_linerange(aes(x = x, ymax = ymax, ymin = ymin), data = line_between_reference_and_undertransposed, linetype = 'dashed') +
  geom_text(aes(x = x+50, y = mean(c(ymax,ymin)), label = round(ymax-ymin, 2)), data = line_between_reference_and_overtransposed) +
  geom_text(aes(x = x+50, y = mean(c(ymax,ymin)), label = round(ymax-ymin, 2)), data = line_between_reference_and_undertransposed)

den <- ggplot(all[all$stat=="fraction",]) + geom_line(aes(x = fragment_length, y = value, color = sample)) + theme_bw() + xlab("Fragment length (bp)") +
  coord_cartesian(xlim=c(0, 500)) +
  guides(color = guide_legend(title = '')) + scale_color_manual(values = colors) +
  ylab("Fraction of reads") + theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.6, 0.6), legend.background = element_blank())

pdf('example_fld_calculation_cumulative.pdf', width = 4, height = 3)
print(cumulative)
dev.off()

pdf('example_fld_calculation_density.pdf', width = 4, height = 3)
print(den)
dev.off()

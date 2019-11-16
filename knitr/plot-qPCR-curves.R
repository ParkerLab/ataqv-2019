#!/usr/bin/env Rscript 
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')

# qPCR data
day_1 <- read.table('2017-08-10.txt', head = T, as.is = T, sep = '\t')
day_2 <- read.table('2017-08-14.txt', head = T, as.is = T, sep = '\t')
day_1 <- tidyr::gather(day_1, key = 'parameters', value = 'norm_rlu', -Cycle)
day_2 <- tidyr::gather(day_2, key = 'parameters', value = 'norm_rlu', -Cycle)

day_1$parameters <- gsub('^X', '', day_1$parameters)
day_2$parameters <- gsub('^X', '', day_2$parameters)
day_1$parameters <- gsub('\\.\\.', '_', day_1$parameters)
day_2$parameters <- gsub('\\.\\.', '_', day_2$parameters)
day_1$parameters <- gsub('\\.$', '', day_1$parameters)
day_2$parameters <- gsub('\\.$', '', day_2$parameters)
day_1$parameters <- gsub('(rep)\\.(\\d)', '\\1\\2', day_1$parameters)
day_2$parameters <- gsub('(rep)\\.(\\d)', '\\1\\2', day_2$parameters)

# the cycle number is 'additional cycles', i.e. on top of the 5 initial cycles
day_1$Cycle <- day_1$Cycle + 5
day_2$Cycle <- day_2$Cycle + 5

# there were four replicates made each day, but we only used the best 3. So remove the omitted libraries
day_1 <- day_1[grep('rep3', day_1$parameters, invert = T),]
day_1$parameters <- gsub('rep4', 'rep3', day_1$parameters)
day_2 <- day_2[grep('rep1', day_2$parameters, invert = T),]
day_2$parameters <- gsub('rep4', 'rep6', day_2$parameters)
day_2$parameters <- gsub('rep3', 'rep5', day_2$parameters)
day_2$parameters <- gsub('rep2', 'rep4', day_2$parameters)

pcr <- bind_rows(day_1, day_2)
pcr$rep <- gsub('.*_rep', 'rep', pcr$parameters)
pcr$tn5 <- gsub('_rep.*', '', pcr$parameters)
pcr$tn5 <- factor(pcr$tn5, levels = c('No.cell.ctrl', '0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X'), ordered = T)

p <- ggplot(pcr) + geom_line(aes(x = Cycle, y = norm_rlu, color = tn5)) +
  scale_color_manual(values=c('green', tn5_colors), labels = c('no cells (control)', '0.2', '0.5', '0.66', '1', '1.5', '2', '5')) +
  guides(color = guide_legend(title = 'Relative [Tn5]')) +
  theme_bw() + ylab('Normalized RLU') + xlab('Total PCR cycles') +
  geom_vline(xintercept = 16, color = 'red', linetype = 'dashed') +
  facet_wrap(~rep)
pdf('tn5_qpcr.pdf', height = 6, width = 8)
print(p)
dev.off()

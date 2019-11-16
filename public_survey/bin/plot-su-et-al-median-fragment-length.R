#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(viridis)
library(ggplot2)

SAMPLE_INFO <- commandArgs(T)[1]
METRICS <- commandArgs(T)[2]

sample_info <- read.table(SAMPLE_INFO, head = T, sep = '\t', as.is=T) %>%
  dplyr::rename(library=experiment) %>%
  dplyr::filter(project=='PRJNA323617.PRJNA341508')
sample_info$condition <- gsub('GSM\\d+: (.*)_rep.*', '\\1', sample_info$title)

flowcells <- bind_rows(lapply(list.files(path = '.',pattern = 'flowcells.txt'),
                              function(f){
                                tmp <- read.table(f, head=F, as.is = F, sep='\t',col.names = c('flowcell'))
                                tmp$run <- gsub('.flowcells.txt', '', basename(f))
                                return(tmp)
                              }))


metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric','value')) %>%
  dplyr::filter(metric=='median_fragment_length') %>%
  dplyr::select(library, value) %>%
  dplyr::rename(median_fragment_length=value) %>%
  dplyr::mutate(median_fragment_length=as.numeric(median_fragment_length))

sample_info <- left_join(sample_info, metrics) %>%
  left_join(flowcells)


KEEP_FLOWCELL <- unique(sample_info$flowcell[sample_info$condition %in% c('E0', 'E1')])
KEEP_LIBRARIES <- unique(sample_info$library[sample_info$flowcell %in% KEEP_FLOWCELL])
sample_info <- sample_info[sample_info$library %in% KEEP_LIBRARIES,]
stopifnot(length(unique(sample_info$run))==length(unique(sample_info$library)))
sample_info <- sample_info[,c('condition', 'flowcell', 'median_fragment_length')]
sample_info$flowcell <- as.numeric(as.factor(sample_info$flowcell))
print(sample_info)

p <- ggplot(sample_info) +
  geom_point(aes(x = condition, y = median_fragment_length, color = as.factor(flowcell))) +
  theme_bw() +
  scale_color_viridis_d() +
  ylab('Median fragment length (bp)') +
  xlab('Condition') +
  guides(color=guide_legend(title='Flowcell'))
pdf('su_et_al_median_fragment_length_vs_batch.pdf', height = 3, width = 5)
p
dev.off()

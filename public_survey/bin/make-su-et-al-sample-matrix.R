#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

SAMPLE_INFO <- commandArgs(T)[1]
METRICS <- commandArgs(T)[2]

sample_info <- read.table(SAMPLE_INFO, head = T, sep = '\t', as.is=T) %>%
  dplyr::rename(library=experiment) %>%
  dplyr::filter(project=='PRJNA323617.PRJNA341508')
sample_info$condition <- gsub('.* (E\\d+)_rep.*', '\\1', sample_info$title)
sample_info <- sample_info[sample_info$condition %in% c('E0', 'E1'),c('run', 'library', 'condition')]

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
write.table(sample_info[,c('median_fragment_length', 'library', 'condition')], 'covariates.with_median_fragment_length.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)
write.table(sample_info[,c('library', 'condition')], 'covariates.no_median_fragment_length.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = T)

#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(optparse)
library(broom)

option_list <- list(
  make_option(c("--sample_info"), action = 'store', type = 'character', help = '[Required] Path to sample info file'),
  make_option(c("--out"), action = 'store', type = 'character', help = '[Required] Name of output file')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

# GLOBALS -- should be shared across plotting scripts
tn5_levels <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

tn5_replicate_colors <- c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
names(tn5_replicate_colors) <- paste0('', 1:6)

tn5_colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')

chromhmm_colors <- c('weak_enhancer' = rgb(255,252,4,maxColorValue=255),
                     'strong_enhancer' = rgb(250,202,0,maxColorValue=255),
                     'transcribed' = rgb(0,176,80,maxColorValue=255),
                     'insulator' = rgb(10,190,254,maxColorValue=255),
                     'weak_promoter' = rgb(255,105,105,maxColorValue=255),
                     'active_promoter' = rgb(255,0,0,maxColorValue=255),
                     'poised_promoter' = rgb(207,11,198,maxColorValue=255),
                     'repressed' = rgb(127,127,127,maxColorValue=255),
                     'low_signal' = rgb(0,0,0,maxColorValue=255))
chromhmm_color_labels <- gsub('_', ' ', names(chromhmm_colors))
names(chromhmm_color_labels) <- names(chromhmm_colors)


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

chromhmm_overlap <- bind_rows(lapply(list.files(full.names = T, pattern = '.counts.bed'), function(f){
  tmp <- read.table(f, head = F, as.is = T)
  colnames(tmp) <- c('read_count', 'chromatin_state')
  tmp$library <- gsub('.counts.bed', '', basename(f))
  return(tmp)}))


counts <- left_join(chromhmm_overlap, sample_info)
counts$tn5 <- log2(as.numeric(gsub('X', '', counts$tn5)))
counts <- counts %>% dplyr::group_by(library) %>% dplyr::mutate(total_reads=sum(read_count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(proportion_of_reads=read_count/total_reads)
merge_reps <- counts %>% dplyr::select(chromatin_state, proportion_of_reads, tn5) %>% dplyr::group_by(tn5, chromatin_state) %>% dplyr::summarize(mean_proportion=mean(proportion_of_reads), median_proportion=median(proportion_of_reads))
# change the order of the chromatin states so that the order in the legend matches the order of the lines on the plot
tmp <- merge_reps[merge_reps$tn5==max(merge_reps$tn5),]
tmp <- tmp[order(tmp$median_proportion),]
chromatin_state_order <- rev(tmp$chromatin_state)
merge_reps$chromatin_state <- factor(merge_reps$chromatin_state, levels = chromatin_state_order, ordered = T)

p <- ggplot(merge_reps) + geom_line(aes(x = tn5, y = median_proportion, color = chromatin_state)) +
  geom_point(aes(x = tn5, y = median_proportion, color = chromatin_state)) +
  scale_color_manual(values=chromhmm_colors, labels = chromhmm_color_labels) +
  theme_bw() + xlab('log2(Relative [Tn5])') + ylab('Median proportion of reads from\nlibrary overlapping chromatin state') +
  guides(color=guide_legend(title = 'Chromatin state')) +
  theme(legend.position = 'left', legend.background = element_rect(color = 'black', fill = NA))
pdf('tn5_chromatin_state_overlap.pdf', width = 5, height = 3)
print(p)
dev.off()

# model proportion_reads_in_state ~ replicate + tn5
# TODO: check this output
models <- bind_rows(lapply(unique(counts$chromatin_state),
                           function(x){
                             for_mod <- counts[counts$chromatin_state==x,]
                             mod <- tidy(lm(formula = proportion_of_reads ~ replicate + tn5, data = for_mod))
                             mod$state <- x
                             return(mod)
                           }))
models <- models %>%
  dplyr::filter(term=='tn5')
models$q <- models$p.value / nrow(models)
print(models)

#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(RColorBrewer)
library(broom)


# need the peak counts
# need peak-TF overlaps
# need list of Tn5 sensitive peaks
args <- commandArgs(T)
SAMPLE_INFO <- args[1]
TN5_SENSITIVITY <- args[2]
MASTER_PEAK_OVERLAP <- args[3]

tn5_sensitivity <- read.table(TN5_SENSITIVITY, head = F, as.is = T, sep = '\t', col.names = c('peak', 'sensitivity'))
master_peak_overlap <- read.table(MASTER_PEAK_OVERLAP, head = F, as.is = T, sep = '\t', col.names = c('chrom', 'start', 'end', 'overlaps', 'tf')) %>%
  dplyr::mutate(peak=glue('{chrom}:{start}:{end}')) %>%
  dplyr::select(-chrom, -end, -start)
master_peak_overlap <- master_peak_overlap[master_peak_overlap$peak %in% tn5_sensitivity$peak,] # exclude peaks where model didn't converge
master_peak_overlap$overlaps[master_peak_overlap$overlaps>1] <- 1

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t') %>%
	dplyr::select(library, tn5) %>% unique()
sample_info$library <- as.character(sample_info$library)

peak_counts <- bind_rows(lapply(list.files(pattern='.counts.bed'), function(f){
			tmp <- read.table(f, head = F, as.is = T, sep = '\t', col.names = c('chrom', 'start', 'end', 'count'))
			tmp$library <- gsub('.counts.bed', '', basename(f))
			return(tmp)
		  }))
peak_counts$peak <- with(peak_counts, paste(chrom, start, end, sep = ':'))
peak_counts <- dplyr::select(peak_counts, peak, library, count)

peak_signal <- peak_counts %>%
  left_join(sample_info) %>%
  dplyr::filter(tn5=='1X') %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(median_count_in_library=median(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(signal=count/median_count_in_library) %>%
  dplyr::group_by(peak) %>%
  dplyr::summarize(median_signal=median(signal))

peak_signal_bins <- peak_signal
peak_signal_bins$bin <- cut_number(peak_signal_bins$median_signal, n = 10)

tmp <- master_peak_overlap %>%
  left_join(tn5_sensitivity) %>%
  left_join(peak_signal_bins) %>%
  dplyr::mutate(tf=gsub('.ENCFF.*', '', tf))
tmp2 <- tmp %>%
  dplyr::group_by(tf, bin, overlaps) %>%
  dplyr::summarize(prob_significant=sum(sensitivity!='insensitive')/n())
tmp2$bin <- as.numeric(tmp2$bin)
tmp2$overlaps <- ifelse(tmp2$overlaps, 'yes', 'no')
p <- ggplot(tmp2) + geom_line(aes(x = bin, y = prob_significant, color = overlaps)) +
  facet_wrap(~tf) + theme_bw() + xlab('ATAC-seq peak read count decile') +
  scale_x_continuous(breaks=1:10) + ylab('Probability peak is Tn5 sensitive') +
  guides(color = guide_legend(title='ATAC-seq peak\noverlaps\nTF ChIP-seq peak'))
pdf('tn5_sensitivity_prob_by_chipseq_overlap_and_peak_signal.pdf', width = 11, height = 11)
p
dev.off()


# model prob_nominally_sign ~ peak_size + overlaps_chipseq_peak
model_results <- function(glm, check_converged = T) {
  # extracts the variable names, coefficients, and p-values for the variables.
  
  results <- data.frame(
    variable = names(coef(glm)),
    coefficient = coef(glm),
    coefficient_se = coef(summary(glm))[,2],
    p = coef(summary(glm))[,4]
  )
  
  if (check_converged & !glm$converged) {
    results$coefficient <- NA
    results$coefficient_se <- NA
    results$p <- NA
  }
  
  return(results)
}


logistic_regression <- bind_rows(lapply(unique(master_peak_overlap$tf), function(tf){
  print(glue('Modeling {tf}'))
  tmp <- master_peak_overlap[master_peak_overlap$tf == tf,] %>%
    left_join(tn5_sensitivity) %>%
    left_join(peak_signal) %>%
    dplyr::mutate(tf=gsub('.ENCFF.*', '', tf),
                  significant=ifelse(sensitivity!='insensitive', 1, 0),
                  overlaps_peak=ifelse(overlaps, 1, 0))
  
  mod <- glm(significant ~ median_signal + overlaps_peak, data = tmp, family=binomial(link="logit"))
  df <- model_results(mod)
  df$tf <- tf
  return(df)
}))

# TODO: quantify
binding_predicts_sensitivity <- logistic_regression %>%
  dplyr::filter(variable=='overlaps_peak') %>%
  dplyr::select(-variable)
binding_predicts_sensitivity$significant <- binding_predicts_sensitivity$p <= (0.05 / nrow(binding_predicts_sensitivity))
binding_predicts_sensitivity$direction <- sign(binding_predicts_sensitivity$coefficient)
binding_predicts_sensitivity[,c('tf', 'significant', 'direction')]
print(binding_predicts_sensitivity[!binding_predicts_sensitivity$significant | binding_predicts_sensitivity$direction<0,c('tf', 'significant', 'direction', 'p')])


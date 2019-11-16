#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

EXAMPLE_PROMOTER_PEAK <- 'chr1:40348854:40349418'
EXAMPLE_ENHANCER_PEAK <- 'chr12:125389290:125389653'

args <- commandArgs(T)

modeling <- bind_rows(lapply(list.files(pattern = 'results.txt', full.names = T), function(f){
  RE <- 'covariates_(.*)_subsample_(.*).results.txt'
  covariates <- gsub(RE, '\\1', basename(f))
  subsampled <- gsub(RE, '\\2', basename(f))
  tmp <- read.table(f, head = T, as.is = T, sep = '\t')
  tmp$covariates <- gsub('___', ',', covariates)
  return(tmp)
}))

# exclude peaks for which the model didn't converge for all covariates
number_covariate_combos <- length(unique(modeling$covariates))
x <- table(modeling$peak)
CONVERGED_PEAKS <- names(x)[x==number_covariate_combos]
results <- modeling[modeling$peak %in% CONVERGED_PEAKS,]

# now adjust the p-values...
results$q <- 1
for(i in unique(results$covariates)) {
  results$q[results$covariates==i] <- p.adjust(results$p[results$covariates==i], method='BH')
}
results$q[is.na(results$q)] <- 1


tn5_sensitive_up_peaks <- results$peak[results$covariate=='replicate' & results$q<=0.05 & results$coefficient>0]
tn5_sensitive_down_peaks <- results$peak[results$covariate=='replicate' & results$q<=0.05 & results$coefficient<0]
tn5_sensitivity <- data.frame(peak=CONVERGED_PEAKS, sensitivity='insensitive')
tn5_sensitivity$peak <- as.character(tn5_sensitivity$peak)
tn5_sensitivity$sensitivity <- as.character(tn5_sensitivity$sensitivity)
tn5_sensitivity$sensitivity[tn5_sensitivity$peak %in% tn5_sensitive_up_peaks] <- 'up'
tn5_sensitivity$sensitivity[tn5_sensitivity$peak %in% tn5_sensitive_down_peaks] <- 'down'
write.table(tn5_sensitivity, file='tn5-sensitivity.txt', append=F, sep='\t', col.names=F, row.names=F, quote=F)

# show the shifts in the Z-score distributions
results$Covariates <- gsub(',', ' + ', results$covariates)
results$Covariates <- gsub('_', ' ', results$Covariates)

number_significant <- results %>% dplyr::group_by(Covariates) %>%
  dplyr::summarize(number_total=n(),
                   number_significant=sum(q<=0.05)) %>%
  dplyr::mutate(percent_significant=round(100 * number_significant / number_total, 1))

modeling_results_converged_in_all <- results


example_peaks <- results[results$peak %in% c(EXAMPLE_ENHANCER_PEAK, EXAMPLE_PROMOTER_PEAK),]
example_peaks$label <- ifelse(example_peaks$peak==EXAMPLE_PROMOTER_PEAK, 'promoter peak from Fig 2', 'enhancer peak from Fig 2')
example_peak_colors <- c('promoter peak from Fig 2' = 'red', 'enhancer peak from Fig 2' = 'blue')

p <- ggplot(results) +
  geom_density(aes(x = coefficient / coefficient_se, fill = Covariates), alpha = 0.3) +
  geom_text(aes(x = 15, y = 0.1, label = glue('{percent_significant}% of peaks\nTn5 sensitive')), data = number_significant, size = 4) +
  theme_bw() +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
  geom_vline(aes(xintercept = coefficient / coefficient_se, color = label), linetype = 'dashed', data = example_peaks) +
  scale_color_manual(values = example_peak_colors, guide=guide_legend(title='')) +
  xlab('Z-statistic for Tn5 coefficient') +
  facet_wrap(~Covariates, ncol = 1)
pdf('tn5_modeling.pdf', width = 7, height = 5)
print(p)
dev.off()


library(optparse)

option_list <- list(
  make_option(c('--bed'), action = 'store', type = 'character', help = '[Required] Path to the list of regions of interest'),
  make_option(c('--resolution'), action = 'store', type = 'numeric', help = '[Required] Size of the bins'),
  make_option(c('--flanking'), action = 'store', type = 'numeric', help = '[Required] Amount of sequence to take on each side (should be a multiple of --resolution)'),
  make_option(c('--prefix'), action = 'store', type = 'character', help = '[Required] Prefix for the output')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

options(scipen = 999)

library(dplyr)
library(tidyr)

# for testing
#opts <- list(
#  'bed' = '/lab/work/porchard/hcr_lcr/results/snp_sensitive_pwm_scanning/work/peak_classification/peak_classification.run.txt',
#  'resolution' = 100,
#  'flanking' = 500
#)

prefix <- function(x, p = opts$prefix) {
  return(paste(p, x, sep = '.'))
}

bed <- read.table(opts$bed, head = F, as.is = T, sep = '\t')
colnames(bed) <- c('chrom', 'start', 'end')
widths <- bed$end - bed$start

if (sum(widths > (opts$flanking * 2)) > 0) {
	warning(paste0(sum(widths > (opts$flanking * 2)), ' regions have width greater than --flanking*2 (= ', opts$flanking * 2, ')'))
}

bed <- bed %>% tidyr::unite(col = bed, sep = ':', chrom, start, end)
bins <- data.frame(start=seq(-1*opts$flanking, opts$flanking-opts$resolution, opts$resolution),
                   end=seq(-1*opts$flanking+opts$resolution, opts$flanking, opts$resolution))
bins <- bins %>% tidyr::unite(col = bin, sep = ' - ', start, end)

features <- expand.grid(bed=bed$bed, bin=bins$bin)
features <- tidyr::separate(features, col = bed, into = c('chrom', 'start', 'end'), sep = ':', convert = T)
features <- tidyr::separate(features, col = bin, into = c('bin_start', 'bin_end'), sep = ' - ', convert = T)
features <- features[order(features$chrom, features$start),]
features <- tidyr::unite(data = features, col=feature, sep = ':', chrom:end, remove = F)
features$center <- round(0.5*(features$end - features$start) + features$start)
features$start <- features$bin_start + features$center
features$end <- features$bin_end + features$center
features <- features[,c('chrom', 'start', 'end', 'feature')]

stopifnot(nrow(features)==nrow(bed) * nrow(bins))
write.table(x = features, file = prefix('features.bed'), quote = F, sep = '\t', row.names = F, col.names = F, append = F)

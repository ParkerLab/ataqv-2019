library(dplyr)
library(tidyr)

# sample information
# ataqv metrics
# choose the non-outliers, with the same cell type as the outliers and some minimum read depth...
OUTLIER <- 'SRX859231'
SAMPLE_INFO <- file.path(Sys.getenv('ATAQV_ROOT'), 'public_survey/sample_info/sample_info.txt')
METRICS <- file.path(Sys.getenv('ATAQV_ROOT'), 'public_survey/work/downstream/results/ataqv-new-metrics/public-single-cell-hg19.metrics.txt')

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t') %>%
  dplyr::filter(is_single_cell) %>%
  dplyr::select(experiment, organism, title) %>%
  unique()
sample_info$cell_type <- gsub('.*singles-(.*?)(-rep\\d+)?-well-.*', '\\1', sample_info$title)
ORGANISM <- sample_info$organism[sample_info$experiment==OUTLIER]
CELL_TYPE <- sample_info$cell_type[sample_info$experiment==OUTLIER]
from_same_species_and_cell_type <- sample_info$experiment[sample_info$organism==ORGANISM & sample_info$cell_type==CELL_TYPE]

metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric', 'value')) %>%
  tidyr::spread(key = metric, value = value) %>%
  dplyr::select(library, hqaa, max_fraction_reads_from_single_autosome) %>%
  dplyr::mutate(hqaa=as.numeric(hqaa),
                max_fraction_reads_from_single_autosome=as.numeric(max_fraction_reads_from_single_autosome))
metrics <- metrics[metrics$library %in% from_same_species_and_cell_type,]
print('Outlier HQAA:')
print(metrics$hqaa[metrics$library==OUTLIER])
lower_bound <- quantile(metrics$max_fraction_reads_from_single_autosome, c(0.3))
upper_bound <- quantile(metrics$max_fraction_reads_from_single_autosome, c(0.7))
nonoutliers <- metrics[metrics$hqaa>=20000 & metrics$max_fraction_reads_from_single_autosome<=upper_bound & metrics$max_fraction_reads_from_single_autosome>=lower_bound,]

write.table(nonoutliers[,c('library')], file = 'nonoutliers.txt', append = F, quote = F, sep = '\t', row.names = F, col.names = F)

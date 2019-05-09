library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)


args <- commandArgs(trailingOnly = T)
METRICS_FILE <- args[1]
SAMPLE_INFO <- args[2]

ROOT <- Sys.getenv('ATAQV_HOME')

source(file.path(ROOT, "src/convert_tn5.R"))

# for testing
#METRICS_FILE <- '/lab/work/porchard/ataqv.final/cluster_density/separate/work/extract_ataqv_metrics/results/cluster_density_tn5_series.metrics.txt'
#SAMPLE_INFO <- '/lab/work/porchard/ataqv.final/cluster_density/separate/sample_info/sample_info.txt'

metrics <- read.table(METRICS_FILE, head = F, as.is = T, sep = '\t')
colnames(metrics) <- c("file", 'metric', 'value')

parse_file_name <- function(f) {
  regex <- "^(\\d+)___(\\d+)\\.json.gz$"
  seqcore_id <- gsub(regex, "\\1", basename(f))
  sequencing_run <- gsub(regex, "\\2", basename(f))
  return(c('seqcore_id' = seqcore_id, 'sequencing_run' = sequencing_run))
}

metrics$seqcore_id <- sapply(metrics$file, function(x){parse_file_name(x)['seqcore_id']})
metrics$sequencing_run <- sapply(metrics$file, function(x){parse_file_name(x)['sequencing_run']})

sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$seqcore_id <- as.character(sample_info$seqcore_id)
sample_info$sequencing_run <- as.character(sample_info$sequencing_run)

metrics <- left_join(metrics, sample_info) %>% dplyr::select(metric, value, seqcore_id, sequencing_run, tn5, replicate)
metrics$sequencing_run <- paste('run', metrics$sequencing_run)

metrics <- metrics %>% tidyr::spread(key = sequencing_run, value = value)
metrics$tn5 <- as.factor.tn5(metrics$tn5)

# set the color scheme
# blues <- seq(1, 0, length.out = 7)
# reds <- rev(seq(1, 0, length.out = 7))
# greens <- rep(0, 7)
# colors <- rep(NA, 7)

# for(i in 1:length(colors)) {
#   colors[i] <- rgb(reds[i], greens[i], blues[i])
# }
# names(colors) <- sort(unique(metrics$tn5))
colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')
names(colors) <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

p <- ggplot(metrics) + geom_point(aes(x = `run 1830`, y = `run 1833`, color = tn5, shape = replicate)) +
  theme_bw() +
  facet_wrap(~metric, scales = 'free') +
  xlab('Low cluster density (sequencing run 1830)') +
  ylab('High cluster density (sequencing run 1833)') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_color_manual(values = colors)
pdf('effect_of_cluster_density_on_ataqv_metrics.pdf', width = 15, height = 10)
print(p)
dev.off()


for(each_metric in unique(metrics$metric)) {
  tmp <- metrics %>% dplyr::filter(metric==each_metric)
  tmp$metric <- dplyr::recode(tmp$metric,
                              hqaa = 'High-quality autosomal reads',
                              percent_mitochondrial = '% mitochondrial reads',
                              median_fragment_length = 'Median fragment length (bp)',
                              total_peaks = 'Number of peaks',
                              hqaa_overlapping_peaks_percent = '% high-quality autosomal\nreads overlapping peaks',
                              tss_enrichment = 'TSS enrichment',
                              short_mononucleosomal_ratio = 'Short : mononucleosomal reads',
                              fragment_length_distance = 'Fragment length distance'
  )
  
  xlabel <- paste(unique(tmp$metric), '\n(low cluster density)')
  ylabel <- paste(unique(tmp$metric), '\n(high cluster density)')
  lim.max <- max(c(tmp[,'run 1830'], tmp[,'run 1833']))
  lim.min <- min(c(tmp[,'run 1830'], tmp[,'run 1833']))
  p <- ggplot(tmp) + geom_point(aes(x = `run 1830`, y = `run 1833`, color = tn5, shape = replicate)) +
    theme_bw() +
    xlab(xlabel) +
    ylab(ylabel) +
    coord_cartesian(xlim=c(lim.min, lim.max),
                    ylim=c(lim.min, lim.max)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    scale_color_manual(values = colors)
  
  f <- paste0('effect_of_cluster_density_on_', each_metric, '.pdf')
  pdf(f, height = 4, width = 5)
  print(p)
  dev.off()
}

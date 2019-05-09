library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(glue)

args <- commandArgs(trailingOnly = T)
METRICS_FILE <- args[1]
SAMPLE_INFO <- args[2]
PREFIX <- args[3]

ROOT <- Sys.getenv('ATAQV_HOME')

source(file.path(ROOT, "src/convert_tn5.R"))

prefix <- function(x) {
  paste(PREFIX, x, sep = '.')
}

# for testing
# ROOT <- '/lab/work/porchard/ataqv.final/cluster_density'
# METRICS_FILE <- file.path(ROOT, 'figures/ataqv_metrics/metrics.txt')
# SAMPLE_INFO <- file.path(ROOT, 'sample_info/sample_info.txt')


sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$library <- sample_info$id
sample_info$tn5 <- as.factor.tn5(sample_info$tn5)

parse_file_name <- function(f) {
  RE <- '^(.*).json.gz'
  if (grepl(RE, basename(f))) {
    lib <- gsub(RE, '\\1', basename(f))
    x <- c('library' = lib)
    return(x)
  } else {
    warning('Cant parse file name!')
    quit(save = 'no')
  }
}


dat <- read.table(METRICS_FILE, head = F, as.is = T, sep = '\t')
colnames(dat) <- c("file", "metric", 'value')

dat$library <- sapply(dat$file, function(x) {parse_file_name(x)['library']})
dat <- left_join(dat, sample_info)


metric_names <- c(hqaa = 'High-quality autosomal reads', 
                  percent_autosomal_duplicate = '% autosomal duplicate', 
                  percent_mitochondrial = '% mitochondrial reads', 
                  median_fragment_length = 'Median fragment length (bp)', 
                  total_peaks = 'Number of peaks', 
                  hqaa_overlapping_peaks_percent = '% high-quality autosomal\nreads overlapping peaks', 
                  tss_enrichment = 'TSS enrichment', 
                  short_mononucleosomal_ratio = 'Short : mononucleosomal reads', 
                  fragment_length_distance = 'Fragment length distance')


# first, plot the cluster density stuff
metrics <- dat %>% dplyr::select(-library, -id, -file)
metrics$sequencing_run <- paste('run', metrics$sequencing_run)

metrics <- metrics %>% tidyr::spread(key = sequencing_run, value = value)
metrics$tn5 <- as.factor.tn5(as.numeric.tn5(metrics$tn5))

# set the color scheme
colors <- c('#0000ff','#2510a3','#1f1050','#000000','#52170b','#a51b0b','#ff0000')
names(colors) <- c('0.2X', '0.5X', '0.66X', '1X', '1.5X', '2X', '5X')

p <- ggplot(metrics) + geom_point(aes(x = `run 1830`, y = `run 1833`, color = tn5, shape = replicate)) +
  theme_bw() +
  facet_wrap(~metric, scales = 'free', labeller = as_labeller(metric_names)) +
  xlab('Low cluster density') +
  ylab('High cluster density') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_color_manual(values = colors)
pdf(prefix('effect_of_cluster_density_on_ataqv_metrics.pdf'), width = 15, height = 10)
print(p)
dev.off()


for(each_metric in unique(metrics$metric)) {
  tmp <- metrics %>% dplyr::filter(metric==each_metric)
  
  xlabel <- paste(metric_names[each_metric], '\n(low cluster density)')
  ylabel <- paste(metric_names[each_metric], '\n(high cluster density)')
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
  
  f <- prefix(glue('effect_of_cluster_density_on_{each_metric}.pdf'))
  pdf(f, height = 4, width = 5)
  print(p)
  dev.off()
}



# next, plot the Tn5 stuff
metrics <- metrics %>% tidyr::gather(key = sequencing_run, value = value, `run 1830`, `run 1833`)

for(sequencing_run in unique(metrics$sequencing_run)) {
  seqrun <- gsub(' ', '_', sequencing_run)
  dat <- metrics[metrics$sequencing_run==sequencing_run,]
  
  p <- ggplot(dat) + geom_point(aes(x = tn5, y = value, color = replicate)) +
    theme_bw() +
    facet_wrap(~metric, scales = 'free', labeller = as_labeller(metric_names)) +
    xlab('Tn5')
  
  f <- prefix(glue('tn5_vs_ataqv_metrics.{seqrun}.pdf'))
  pdf(f, height = 7, width = 10)
  print(p)
  dev.off()
  
  for(each_metric in unique(dat$metric)) {
    tmp <- dat %>% dplyr::filter(metric==each_metric)
    
    p <- ggplot(tmp) + geom_point(aes(x = tn5, y = value, color = replicate)) +
      theme_bw() +
      # scale_color_manual(values=colors) +
      theme(axis.text=element_text(size=13))+
      theme(axis.title.y=element_text(size=13))+
      xlab("Tn5") +
      ylab(metric_names[each_metric])
    
    f <- prefix(glue('tn5_vs_{each_metric}.{seqrun}.pdf'))
    pdf(f, height = 5, width = 6)
    print(p)
    dev.off()
    
    # p <- p + geom_text_repel(aes(x = tn5, y = value, label = seqcore_id))
    # f <- prefix(paste0('effect_of_tn5_concentration_on_', each_metric, '.seqcore_id_labels.pdf'))
    # pdf(f, height = 5, width = 6)
    # print(p)
    # dev.off()
  }
}

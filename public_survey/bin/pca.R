#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--counts-dir"), dest = 'counts_dir', action = 'store', type = 'character', help = '[Required] Path to the directory containing *.counts.bed files'),
  make_option(c("--scale"), action = 'store_true', default = F, type = 'logical', help = '[Optional] scale before PCA'),
  make_option(c("--log2"), action = 'store_true', default = F, type = 'logical', help = '[Optional] Take log2 of the (RPKM matrix + pseudocount)'),
  make_option(c("--pseudocount"), action = 'store', default = 1, type = 'numeric', help = '[Optional] Pseusocount to add in the case that log2 of the RPKM matrix is being used'),
  make_option(c("--prefix"), action = 'store', type = 'character', help = '[Required] Prefix for the ouput files')
)

description <- 'Perform PCA on the master peak counts and save the resulting prcomp object to a file.'

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = description)
opts <- parse_args(option_parser)

PREFIX <- opts$prefix

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(glue))

source("~/src/counts_to_rpkm.R")


# load in the counts
counts_files <- list.files(opts$counts_dir, pattern = '\\.bed$', full.names = T)

parse_file_name <- function(f) {
  regex <- '^(.*)\\.bed$'
  id <- gsub(regex, "\\1", basename(f))
  x <- c('id' = id)
  return(x)
}

load_counts_file <- function(f) {
  x <- read.table(f, head = F, as.is = T, sep = '\t')
  colnames(x) <- c('chrom', 'start', 'end', 'count')
  x$id <- parse_file_name(f)['id']
  return(x)
}

counts <- lapply(counts_files, load_counts_file)
counts <- bind_rows(counts)
counts <- counts %>% tidyr::spread(key=id, value=count)

region_lengths <- counts$end - counts$start
rownames(counts) <- with(counts, paste(chrom, start, end, sep = ':'))
counts <- counts %>% dplyr::select(-chrom, -start, -end)

# convert to RPKM
rpkm <- apply(counts, 2, function(x){counts_to_rpkm(x, lengths = region_lengths)})

if (opts$log2) {
  rpkm <- log2(rpkm + opts$pseudocount)
}

# perform PCA
pca <- prcomp(t(rpkm), center = T, scale. = opts$scale)

# write the variances
variances <- pca$sdev^2
variances <- variances / sum(variances)
variances <- data.frame('PC'=1:length(variances), proportion_of_variance=variances)
write.table(variances, file = glue('{PREFIX}variances.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# write the scores
scores <- as.data.frame(pca$x)
scores$library <- rownames(scores)
write.table(scores, file = glue('{PREFIX}scores.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)

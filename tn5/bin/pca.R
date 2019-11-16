#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--counts-dir"), dest = 'counts_dir', action = 'store', type = 'character', help = '[Required] Path to the directory containing *.counts.bed files'),
  make_option(c("--scale"), action = 'store_true', default = F, type = 'logical', help = '[Required] scale before PCA'),
  make_option(c("--out"), action = 'store', type = 'character', help = '[Required] Name of the .Rda file that will contain the pca object')
)

description <- 'Perform PCA on the master peak counts and save the resulting prcomp object to a file.'

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = description)
opts <- parse_args(option_parser)

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

counts_to_rpkm <- function(counts_vector, lengths) {
  # lengths indicates the lengths (in bp) of each of the corresponding regions
  stopifnot(length(counts_vector) == length(lengths))
  
  counts_sum <- sum(counts_vector)
  
  # to reduce the probability of integer overflow,
  # enforce a certain order of operations
  rpkm <- counts_vector * (((10^9) / (lengths)) / counts_sum)
  
  return(rpkm)
}

# load in the counts
counts_files <- list.files(opts$counts_dir, pattern = '\\.counts\\.bed$', full.names = T)

parse_file_name <- function(f) {
  regex <- '^(.*)\\.counts.bed$'
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

# perform PCA
pca <- prcomp(t(rpkm), center = T, scale. = opts$scale)

save(pca, file = opts$out)

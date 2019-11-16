#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("--counts"), action = 'store', type = 'character', help = '[Required] Path to the directory of peak read counts'),
  make_option(c("--out"), action = 'store', type = 'character', help = '[Required] File to which results will be written'),
  make_option(c("--covariates"), action = 'store', type = 'character', help = '[Optional] List of covariates to include. Either ataqv metrics, or PCs (example format: fragment_length_distance,PC1)')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)


# suppressPackageStartupMessages(library("MASS")) # for GLM
suppressPackageStartupMessages(library("DESeq2")) # for GLM
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library("dplyr"))

set.seed(77689)

# load in the peak read counts for each experiment
count_files <- list.files(opts$counts, pattern = 'counts.bed$', full.names = T)

parse_filename <- function(f) {
  info <- gsub(".counts.bed", "", basename(f))
  names(info) <- c('library')
  return(info)
}

load_count_file <- function(f) {
  tmp <- read.table(f, header = F, as.is = T, sep = '\t')
  colnames(tmp) <- c('chrom', 'start', 'end', 'count')
  tmp$library <- parse_filename(f)
  return(tmp)
}

counts <- bind_rows(lapply(count_files, load_count_file))
counts <- counts %>% tidyr::spread(key = library, value = count)
stopifnot(all(!is.na(counts)))

rownames(counts) <- with(counts, paste(chrom, start, end, sep = ':'))
counts <- counts %>% dplyr::select(-chrom, -start, -end)

# set up for DESeq2
exp_info <- read.table(opts$covariates, head = T, as.is = T, sep = '\t')
rownames(exp_info) <- exp_info$library
exp_info <- dplyr::select(exp_info, -library)

# put variable of interest at the end...
VARIABLE_OF_INTEREST <- 'condition'
exp_info_cols <- colnames(exp_info)
exp_info_cols <- c(exp_info_cols[exp_info_cols!=VARIABLE_OF_INTEREST], VARIABLE_OF_INTEREST)
exp_info <- exp_info[,exp_info_cols,drop=F]

# only keep the libraries that we will use...
exp_info <- exp_info[rownames(exp_info) %in% colnames(counts),,drop=F]

des <- formula(paste0('~ ', paste(colnames(exp_info), collapse = ' + ')))
print('Design:')
print(des)
dds <- DESeqDataSetFromMatrix(countData = counts[,rownames(exp_info)], colData=exp_info, design=des)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, quiet=T)
print('Finished dispersion estimates')
dds <- nbinomWaldTest(dds, betaPrior = T)
print('Wald test')
res <- as.data.frame(results(dds, alpha = 0.05))
res$peak <- rownames(res)

# write out the results
# downstream processing will be done with python
write.table(res, opts$out, append = F, quote = F, col.names = T, sep = '\t', row.names = F)

#!/usr/bin/env Rscript
# Build a NB GLM for Tn5-sensitive peaks

source(file.path(Sys.getenv('ATAQV_HOME'), 'src/convert_tn5.R'))
suppressPackageStartupMessages(library('optparse'))

option_list <- list(
  make_option(c('--sample_info'), action = 'store', type = 'character', help = '[Required] Path to sample_info.txt'),
  make_option(c('--counts'), action = 'store', type = 'character', help = '[Required] Path to the directory of peak read counts'),
  make_option(c('--prefix'), action = 'store', type = 'character', help = '[Required] Prefix to use for the outfiles'),
  make_option(c('--covariates'), action = 'store', type = 'character', help = '[Optional] List of covariates to include (ataqv metrics) (example format: fragment_length_distance,tss_enrichment)'),
  make_option(c('--covariates_file'), action = 'store', type = 'character', help = '[Optional] Path to the file of extracted ataqv covariates (required if ataqv metrics are used as covariates)'),
  make_option(c('--tn5_scale'), action = 'store', type = 'character', help = '[Required] "linear" or "log2"'),
  make_option(c('--size_factors'), action = 'store', type = 'character', default = 'hqaa', help = '[Optional] Use this covariate as the size factors (otherwise, use the # reads in peaks as calculated from the counts files as the size factors)')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts <- parse_args(option_parser)

counts_to_rpkm <- function(counts_vector, lengths) {
  # lengths indicates the lengths (in bp) of each of the corresponding regions
  stopifnot(length(counts_vector) == length(lengths))
  
  counts_sum <- sum(counts_vector)
  
  # to reduce the probability of integer overflow,
  # enforce a certain order of operations
  rpkm <- counts_vector * (((10^9) / (lengths)) / counts_sum)
  
  return(rpkm)
}


suppressPackageStartupMessages(library('MASS'))
suppressPackageStartupMessages(library('tidyr'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('dplyr'))

set.seed(77689)

prefix <- function(x, p = opts$prefix) {
  paste(p, x, sep='.')
}

covariates <- strsplit(opts$covariates, ',')[[1]]

sample_info <- read.table(opts$sample_info, head = T, as.is = T, sep = '\t') %>% dplyr::select(library, tn5, replicate) %>% unique() %>%
  dplyr::mutate(library=as.character(library))


if (!is.null(opts$covariates_file)) {
	covariate_values <- read.table(opts$covariates_file, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric', 'value')) %>%
	  dplyr::mutate(library=as.character(library)) %>%
	  tidyr::spread(key = metric, value = value)
	sample_info <- left_join(sample_info, covariate_values)
}

if (opts$tn5_scale == 'log2') {
  sample_info$tn5 <- log2(as.numeric.tn5(sample_info$tn5))
}

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
  tmp$library <- parse_filename(f)['library']
  return(tmp)
}

counts <- bind_rows(lapply(count_files, load_count_file)) %>% tidyr::spread(key = library, value = count)

rownames(counts) <- with(counts, paste(chrom, start, end, sep = ':'))
counts <- counts %>% dplyr::select(-chrom, -start, -end)

size_factors <- data.frame(library = colnames(counts),
                           size = colSums(counts))

if(!is.null(opts$size_factors) & opts$size_factors %in% colnames(sample_info)) {
  size_factors <- sample_info[,c('library', opts$size_factors)]
  colnames(size_factors) <- c('library', 'size')
}

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


# set it up to run over the entire data.frame
process_row <- function(row_number) {
  for_model <- tidyr::gather_(counts[row_number,], key = "library", value = "count", colnames(counts), convert = F)
  for_model <- left_join(for_model, size_factors)
  peak_of_interest <- rownames(counts)[row_number]
  for_model <- left_join(for_model, unique(sample_info))
  
  # model formula
  form <- paste("count ~", paste(covariates, collapse = ' + '), '+ tn5')
  mod <- NULL
  form <- paste(form, '+ offset(log(size))')
  
  tryCatch(
    {
      mod <- glm.nb(formula = as.formula(form), data = for_model)
    }, warning = function(w) {
      # skip
      m <- paste("Warning for peak ", peak_of_interest, ": ", w, "; skipping", sep = '')
      print(m)
      return("")
    }, error = function(e) {
      # Skip
      m <- paste("Error for peak ", peak_of_interest, ": ", e, "; skipping", sep = '')
      print(m)
      return("")
    }
  )

  # extract results

  if(is.null(mod)) {
    return("")
  } else {
    results <- model_results(mod)
    results$variable <- gsub("\\(Intercept\\)", "intercept", results$variable)
    results$peak <- peak_of_interest
    
    return(results)
  }
}


x <- list()
for(i in 1:nrow(counts)) {
	tmp <- process_row(i)
	if (is.data.frame(tmp)) {
		x[[length(x) + 1]] <- tmp
	}
}


all_results <- bind_rows(x)
rownames(all_results) <- as.character(1:nrow(all_results))

simple_results <- all_results %>% dplyr::filter(variable=='tn5') %>%
  dplyr::select(coefficient, coefficient_se, p, peak)

f <- paste(opts$prefix, 'results.txt', sep = '.')
write.table(x = simple_results, file = f, append = F, quote = F, sep = "\t", row.names = F, col.names = T)

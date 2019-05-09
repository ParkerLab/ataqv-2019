library(glue)

args <- commandArgs(T)
# args <- c('/lab/work/porchard/ataqv.final/tn5/merged/work/tn5_sensitive_peaks/modeling/covariates_replicate_subsample_false.results.simple.txt')
modeling_results <- args[1]

x <- read.table(modeling_results, head = T, as.is = T, sep = '\t')
x$q <- p.adjust(x$p, method = 'BH')
number_significant_down <- nrow(x[x$q<=0.05 & x$coefficient<0,])
number_significant_up <- nrow(x[x$q<=0.05 & x$coefficient>0,])
number_significant <- number_significant_up + number_significant_down
number_converged <- nrow(x)

print(glue('Number converged: {number_converged}'))
print(glue('Number significant (FDR 5%): {number_significant}'))
print(glue('Number significant up (FDR 5%): {number_significant_up}'))

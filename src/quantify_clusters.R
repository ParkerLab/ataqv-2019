library(glue)

PF_CLUSTER_FILES <- commandArgs(T)

for(f in PF_CLUSTER_FILES) {
  sequencing_run <- gsub('run_(\\d+).txt', '\\1', basename(f))
  clusters <- sum(read.table(f, head = T, as.is = T, sep = '\t')$PF.Clusters)
  print(glue('There are {clusters} PF clusters in Run {sequencing_run}'))
}
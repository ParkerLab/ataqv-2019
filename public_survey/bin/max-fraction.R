#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(viridis)

args <- commandArgs(T)
SAMPLE_INFO <- args[1]
METRICS <- args[2]
CHROM_COUNTS <- args[3]
PREFIX <- args[4]
BULK_OR_SINGLE_CELL <- args[5]
PROJECT_COLORS <- args[6]

project_colors <- read.table(PROJECT_COLORS, head = T, sep = '\t', as.is = T, comment.char = '')
COLORS <- project_colors$color
names(COLORS) <- project_colors$project

chrom_counts <- read.table(CHROM_COUNTS, head = F, as.is = T, sep = '\t', col.names = c('library', 'chrom', 'reads'))
chrom_counts <- chrom_counts[grepl('^chr\\d+$', chrom_counts$chrom),]
sorted_chromosomes <- paste0('chr', sort(as.numeric(unique(gsub('chr', '', chrom_counts$chrom)))))
chrom_counts$chrom <- factor(chrom_counts$chrom, levels=sorted_chromosomes, ordered = T)
favored_chromosome <- chrom_counts %>%
  dplyr::group_by(library) %>%
  dplyr::mutate(max_count=max(reads)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(reads==max_count) %>%
  dplyr::select(library, chrom) %>%
  dplyr::rename(favored_chrom=chrom)

metrics <- read.table(METRICS, head = F, as.is = T, sep = '\t', col.names = c('library', 'metric', 'value')) %>%
  dplyr::mutate(value=as.numeric(value))

# assign cell type to all public libraries...
sample_info <- read.table(SAMPLE_INFO, head = T, as.is = T, sep = '\t')
sample_info$library <- sample_info$experiment
sample_info$cell_type <- NA
sample_info$cell_type[sample_info$project=='PRJNA259243'] <- 'T cells'
sample_info$cell_type[grepl('H1ESC', sample_info$title, ignore.case = T)] <- 'H1 ESC'
sample_info$cell_type[grepl(' T cell', sample_info$title, ignore.case = T)] <- 'T cells'
sample_info$cell_type[grepl(' B cell', sample_info$title, ignore.case = T)] <- 'B cells'
sample_info$cell_type[grepl('GM12878', sample_info$title)] <- 'GM12878'
sample_info$cell_type[grepl('-GM-', sample_info$title)] <- 'GM12878'
sample_info$cell_type[grepl('K562', sample_info$title, ignore.case = T)] <- 'K562'
sample_info$cell_type[grepl('-TF1-', sample_info$title)] <- 'TF1'
sample_info$cell_type[grepl('HL60', sample_info$title)] <- 'HL60'
sample_info$cell_type[grepl('-BJ-', sample_info$title)] <- 'BJ'
sample_info$cell_type[grepl('-EML-', sample_info$title)] <- 'EML'
sample_info$cell_type[grepl('mESC', sample_info$title, ignore.case = T)] <- 'mESC'
sample_info$cell_type[grepl('NPCs', sample_info$title, ignore.case = T)] <- 'NPCs'
sample_info$cell_type[sample_info$experiment=='SRX766230'] <- 'Primary proliferating myoblasts'
sample_info$cell_type[sample_info$project=='PRJNA279685'] <- 'primary human neonatal keratinocytes'
sample_info$cell_type[sample_info$project=='PRJNA290920'] <- 'B cells'
sample_info$cell_type[sample_info$project=='PRJNA261246'] <- 'MLL-AF9 leukemic cells'
sample_info$cell_type[sample_info$project=='PRJNA296371'] <- 'CD34+ HSPC'
sample_info$cell_type[sample_info$project=='PRJNA306754' & grepl('Alpha', sample_info$title)] <- 'Pancreatic alpha cells'
sample_info$cell_type[sample_info$project=='PRJNA306754' & grepl('Beta', sample_info$title)] <- 'Pancreatic beta cells'
sample_info$cell_type[sample_info$project=='PRJNA306754' & grepl('Acinar', sample_info$title)] <- 'Pancreatic acinar cells'
sample_info$cell_type[sample_info$project=='PRJNA323617.PRJNA341508'] <- 'dentate granule neurons'
sample_info$cell_type[sample_info$project=='PRJNA325538' & grepl('_rod_', sample_info$title)] <- 'Rods'
sample_info$cell_type[sample_info$project=='PRJNA325538' & grepl('Cone', sample_info$title)] <- 'Cones'
sample_info$cell_type[sample_info$project=='PRJNA325538' & grepl('Photoreceptor', sample_info$title)] <- 'Photoreceptor'
sample_info$cell_type[sample_info$project=='PRJNA319573' & grepl('H7 human embryonic', sample_info$title)] <- 'H7 hESC-derived'
sample_info$cell_type[sample_info$project=='PRJNA326967'] <- 'Fetal skin fibroblasts'
sample_info$cell_type[sample_info$project=='PRJNA330723' & grepl('ESC', sample_info$title)] <- 'mESC'
sample_info$cell_type[sample_info$project=='PRJNA330723' & grepl('NPC', sample_info$title)] <- 'NPCs'
sample_info$cell_type[sample_info$project=='PRJNA345167'] <- 'Mouse cortex'
sample_info$cell_type[sample_info$project=='PRJNA349462'] <- 'T cells'
sample_info$cell_type[sample_info$project=='PRJNA352791' & grepl(' IEL', sample_info$title)] <- 'lymphocytes'
sample_info$cell_type[sample_info$project=='PRJNA358685'] <- 'K562'
sample_info$cell_type[sample_info$project=='PRJNA369195' & grepl(' DE', sample_info$title)] <- 'Endoderm'
sample_info$cell_type[sample_info$project=='PRJNA369195' & grepl(' ESC', sample_info$title)] <- 'mESC'
sample_info$cell_type[sample_info$project=='PRJNA371831'] <- 'Trophoblast stem cell-derived'
sample_info$cell_type[sample_info$project=='PRJNA380816'] <- 'Pre-B cells'
sample_info$cell_type[sample_info$project=='PRJNA380283' & grepl('Primary Frozen Mouse', sample_info$title)] <- gsub('.*Frozen Mouse (.*) using.*', '\\1', sample_info$title[sample_info$project=='PRJNA380283' & grepl('Primary Frozen Mouse', sample_info$title)])
sample_info$cell_type[sample_info$project=='PRJNA380283' & grepl('Primary Frozen Human', sample_info$title)] <- gsub('.*Frozen Human (.*) using.*', '\\1', sample_info$title[sample_info$project=='PRJNA380283' & grepl('Primary Frozen Human', sample_info$title)])
sample_info$cell_type[grepl('Cerebellum', sample_info$title)] <- 'Cerebellum'
sample_info$cell_type[grepl('Corpus Callosum', sample_info$title)] <- 'Corpus Callosum'
sample_info$cell_type[grepl('Frontal Gyrus', sample_info$title)] <- 'Frontal Gyrus'
sample_info$cell_type[grepl('Hippocampus', sample_info$title)] <- 'Hippocampus'
sample_info$cell_type[grepl('Caudate Nucleus', sample_info$title)] <- 'Caudate Nucleus'
sample_info$cell_type[grepl(' T24 ', sample_info$title)] <- 'T24 cells (human bladder cancer)'
sample_info$cell_type[grepl(' Jurkat ', sample_info$title)] <- 'Jurkat cells'
sample_info$cell_type[grepl(' HT1080 ', sample_info$title)] <- 'HT1080'
sample_info$cell_type[grepl(' 293T ', sample_info$title)] <- '293T'
sample_info$cell_type[grepl(' Keratinocytes ', sample_info$title)] <- 'Keratinocytes'
sample_info$cell_type[grepl('beta TC6 and alpha TC1', sample_info$title)] <- 'Pancreatic alpha and beta cell lines'
sample_info$cell_type[sample_info$project=='PRJNA394116' & grepl('MEF', sample_info$title)] <- 'Embryonic fibroblasts'
sample_info$cell_type[grepl('Neurons', sample_info$title)] <- 'Neurons'
sample_info$cell_type[grepl('HCT116', sample_info$title)] <- 'HCT116'
sample_info$cell_type[grepl('TOV21G', sample_info$title)] <- 'TOV21G'
sample_info$cell_type[sample_info$project=='PRJNA412907'] <- 'AVM cells'
sample_info$cell_type[grepl('CD4+', sample_info$title)] <- 'T cells'
# head(sample_info[is.na(sample_info$cell_type),])


if (BULK_OR_SINGLE_CELL=='bulk') {
  metrics <- left_join(metrics, sample_info %>% dplyr::select(library, organism, project, is_single_cell, cell_type) %>% unique()) %>%
    tidyr::spread(key = metric, value = value) %>%
    left_join(favored_chromosome)

  p <- ggplot(metrics) +
    geom_point(aes(x = cell_type, y = max_fraction_reads_from_single_autosome, color=project), alpha=0.3) +
    theme_bw() +
    xlab('Cell type') +
    ylab('Max fraction reads from single autosome') +
    #scale_color_viridis_d(option='plasma') +
    scale_color_manual(values=COLORS) +
    guides(color=guide_legend(title='Project')) +
    coord_flip()
  pdf(glue('{PREFIX}hqaa-vs-max-fraction-group-by-cell-type.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  p <- ggplot(metrics) +
    geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome, color=favored_chrom), alpha=0.5) +
    facet_wrap(cell_type~project) +
    scale_x_log10() +
    theme_bw() +
    xlab('Reads after filtering') +
    ylab('Max fraction reads from single autosome') +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Chromosome'))
  pdf(glue('{PREFIX}hqaa-vs-max-fraction.color-by-chrom.pdf'), height = 10, width = 10)
  print(p)
  dev.off()
} else {
  metrics <- left_join(metrics, sample_info %>% dplyr::select(library, organism, project, is_single_cell, cell_type) %>% unique()) %>%
    tidyr::spread(key = metric, value = value) %>%
    left_join(favored_chromosome)
  p <- ggplot(metrics) +
    geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome, color=cell_type), alpha=0.3) +
    facet_wrap(~project) +
    scale_x_log10() +
    theme_bw() +
    xlab('Reads after filtering') +
    ylab('Max fraction reads from single autosome') +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Cell type'))
  pdf(glue('{PREFIX}hqaa-vs-max-fraction.color-by-cell-type.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  p <- ggplot(metrics) +
    geom_jitter(aes(x = cell_type, y = max_fraction_reads_from_single_autosome, color=favored_chrom), alpha=0.3, height=0, width=0.2) +
    theme_bw() +
    xlab('Cell type') +
    ylab('Max fraction reads from single autosome') +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Chromosome'))
  pdf(glue('{PREFIX}hqaa-vs-max-fraction.color-by-cell-type.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  
  p <- ggplot(metrics) +
    geom_point(aes(x = hqaa, y = max_fraction_reads_from_single_autosome, color=favored_chrom), alpha=0.5) +
    facet_wrap(~project) +
    scale_x_log10() +
    theme_bw() +
    xlab('Reads after filtering') +
    ylab('Max fraction reads from single autosome') +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Chromosome'))
  pdf(glue('{PREFIX}hqaa-vs-max-fraction.color-by-chrom.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  p <- ggplot(metrics) +
    geom_point(aes(x = tss_enrichment, y = max_fraction_reads_from_single_autosome, color=favored_chrom),alpha=0.3) +
    facet_wrap(~project) +
    scale_x_log10() +
    theme_bw() +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Chromosome'))
  pdf(glue('{PREFIX}tss-enrichment-vs-max-fraction.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  p <- ggplot(metrics) +
    geom_point(aes(x = median_fragment_length, y = max_fraction_reads_from_single_autosome, color=favored_chrom), alpha=0.3) +
    facet_wrap(~project) +
    theme_bw() +
    scale_color_viridis_d() +
    guides(color=guide_legend(title='Chromosome'))
  pdf(glue('{PREFIX}mfl-vs-max-fraction.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  outliers <- metrics[metrics$max_fraction_reads_from_single_autosome>=0.15 & metrics$hqaa >= 10000,]
  write.table(outliers, file=glue('{PREFIX}outliers.txt'), append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}


library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(ggrepel)
library(RColorBrewer)
library(broom)

CSV <- commandArgs(T)

nuclear_isolation <- read.csv(CSV, col.names = c('replicate', 'initial_cell_count', 'nuclei_concentration_1', 'nuclei_concentration_2', 'average_nuclei_concentration', 'average_nuclei_count', 'efficiency'))
for(i in colnames(nuclear_isolation)) {
  nuclear_isolation[,i] <- as.numeric(gsub(',', '', nuclear_isolation[,i]))
}
nuclear_isolation <- dplyr::select(nuclear_isolation, replicate, initial_cell_count, average_nuclei_count, efficiency)

# Nuclear isolation efficiency
p <- ggplot(nuclear_isolation) + 
  geom_bar(aes(x = replicate, y = efficiency), stat = 'identity') + 
  theme_bw() +
  ylab('Efficiency of nuclear isolation') + xlab('Replicate') +
  scale_x_continuous(breaks = sort(unique(nuclear_isolation$replicate))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 0.6))
pdf('nuclear_isolation_efficiency.pdf', height = 3, width = 4)
print(p)
dev.off()


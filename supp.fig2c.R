source("permutation_test.R")
library(ggplot2)

set.seed(1234)
data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

k562_rese_granges = readPeakFile(paste(data_path, 'K562_all.txt', sep=''))
k562_faire = readPeakFile(paste(data_path, 'k562_faire.bed.gz', sep=''))
k562_rese_granges$faire = !is.na(findOverlaps(k562_rese_granges, k562_faire, select = 'first'))

k562_rese_granges = subset(k562_rese_granges, faire)

marks = c('H2az', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k79me2', 'H3k9ac', 'H3k9me1', 'H3k9me3', 'H4k20me1')

k562_permutation = data.frame(row.names = marks)
k562_backgrounds = data.frame(matrix(nrow = 20000, ncol = length(marks)))
colnames(k562_backgrounds) = marks

for (mark in marks){
  file = paste(data_path, '/ucsc_k562/', 'wgEncodeBroadHistoneK562', mark, 'StdPk.broadPeak.gz', sep='')
  results = histone_permutation_test(subset(k562_rese_granges, k562_rese_granges$FDR < 0.01), subset(k562_rese_granges, k562_rese_granges$FDR >= 0.01), file)
  k562_permutation[mark, 'p_val'] = results$p_val
  k562_permutation[mark, 'cover_ratio'] = results$cover_ratio
  k562_permutation[mark, 'fold_change'] = results$fold_change
  k562_backgrounds[, mark] = results$backround
}

k562_permutation$q_val = p.adjust(k562_permutation$p_val, method = 'BH')
k562_permutation$name = c('H2A.Z', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me1', 'H3K9me3', 'H4K20me1')

plt = ggplot(k562_permutation, aes(x = name, y = fold_change))+
  geom_point(aes(size = cover_ratio, color = q_val < 0.05))+
  theme_bw() + labs(x = '', color = NULL, size = 'Cover ratio', y = 'Fold enrichment')+
  guides(x =  guide_axis(angle = 45))+
  theme(text = element_text(size = 14))+
  scale_color_manual(labels = c('FDR ≥ 0.05', 'FDR < 0.05'), values = c("black", "red3"))+
  scale_size_continuous(limits = c(0.02, 1))

ggsave('supp.fig2c.png', plot = plt, width = 5, height = 4, units = 'in')

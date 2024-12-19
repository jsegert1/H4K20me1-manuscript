source("permutation_test.R")
library(ggplot2)

set.seed(1234)
data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

hepg2_rese_granges = readPeakFile(paste(data_path, 'HepG2_all.txt', sep=''))
hepg2_faire = readPeakFile(paste(data_path, 'hepg2_faire.bed.gz', sep=''))
hepg2_rese_granges$faire = !is.na(findOverlaps(hepg2_rese_granges, hepg2_faire, select = 'first'))

hepg2_rese_granges = subset(hepg2_rese_granges, faire)

marks = c('H2az', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k79me2', 'H3k9ac', 'H3k9me3', 'H4k20me1')

hepg2_permutation = data.frame(row.names = marks)
hepg2_backgrounds = data.frame(matrix(nrow = 20000, ncol = length(marks)))
colnames(hepg2_backgrounds) = marks

for (mark in marks){
  file = paste(data_path, '/ucsc_hepg2/', 'wgEncodeBroadHistoneHepg2', mark, 'StdPk.broadPeak.gz', sep='')
  results = histone_permutation_test(fg_elements = subset(hepg2_rese_granges, hepg2_rese_granges$FDR < 0.01), bg_elements = subset(hepg2_rese_granges, hepg2_rese_granges$FDR >= 0.01), annotationFile = file)
  hepg2_permutation[mark, 'p_val'] = results$p_val
  hepg2_permutation[mark, 'cover_ratio'] = results$cover_ratio
  hepg2_permutation[mark, 'fold_change'] = results$fold_change
  hepg2_backgrounds[, mark] = results$backround
}

hepg2_permutation$q_val = p.adjust(hepg2_permutation$p_val, method = 'BH')
hepg2_permutation$name = c('H2A.Z', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1')

plt = ggplot(hepg2_permutation, aes(x = name, y = fold_change))+
  geom_point(aes(size = cover_ratio, color = q_val < 0.05))+
  theme_bw() + labs(x = '', color = NULL, size = 'Cover ratio', y = 'Fold enrichment')+
  guides(x =  guide_axis(angle = 45))+
  theme(text = element_text(size = 14))+
  scale_color_manual(labels = c('FDR ≥ 0.05', 'FDR < 0.05'), values = c("black", "red3"))+
  scale_size_continuous(limits = c(0.02, 1))

ggsave('supp.fig2d.png', plot = plt, width = 5, height = 4, units = 'in')

source("permutation_test.R")
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
source('auroc.R')
library(DescTools)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

k562_rese_granges = readPeakFile(paste(data_path, 'K562_all.txt', sep=''))

k562_rese_anno = annotatePeak(k562_rese_granges, TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)

k562_rese_granges$distanceToTSS = abs(data.frame(k562_rese_anno)$distanceToTSS)

k562_rese_granges = subset(k562_rese_granges, distanceToTSS > 5000)

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
  scale_color_manual(labels = c('FDR >= 0.05', 'FDR < 0.05'), values = c("black", "red3"))

ggsave('supp.fig3d.png', plot = plt, width = 5, height = 4, units = 'in')




roc_summary = data.frame()

for (mark in marks){
  file = paste(data_path, '/ucsc_k562/', 'wgEncodeBroadHistoneK562', mark, 'StdPk.broadPeak.gz', sep='')
  results = compute_roc(k562_rese_granges, file)
  results$mark = mark
  roc_summary = rbind(roc_summary, results)
}


k562_aurocs = data.frame(mark = marks, auroc = 0)

for (mark in marks){
  mark_roc = roc_summary[which(roc_summary$mark == mark), ]
  auroc = AUC(mark_roc$FPR, mark_roc$TPR)
  k562_aurocs[which(k562_aurocs$mark == mark), 'auroc'] = auroc
}

k562_aurocs$names = c('H2A.Z', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me1', 'H3K9me3', 'H4K20me1')

k562_aurocs$disp = paste(k562_aurocs$names, round(k562_aurocs$auroc, 3), sep = ' = ')
#k562_aurocs = k562_aurocs[rev(order(k562_aurocs$auroc)),]

plt = ggplot(roc_summary, aes(x = FPR, y = TPR, color = mark)) +
  scale_color_discrete(labels = k562_aurocs$disp)+
  geom_line(linewidth = 1.1)+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  theme(text = element_text(size = 14))+
  coord_cartesian(ylim = c(0, 1), xlim = c(0,1), expand = FALSE)+
  theme(legend.title=element_blank())+
  labs(y= 'True Positive Rate', x = 'False Positive Rate' )

ggsave('supp.fig6e.png', plot = plt, width = 6, height = 4, units = 'in')
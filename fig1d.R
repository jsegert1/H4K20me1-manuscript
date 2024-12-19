source('auroc.R')
library(ggplot2)
library(DescTools)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'
k562_rese = readPeakFile(paste(data_path, 'K562_all.txt', sep=''))
marks = c('H2az', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k79me2', 'H3k9ac', 'H3k9me1', 'H3k9me3', 'H4k20me1')

roc_summary = data.frame()

for (mark in marks){
  file = paste(data_path, '/ucsc_k562/', 'wgEncodeBroadHistoneK562', mark, 'StdPk.broadPeak.gz', sep='')
  results = compute_roc(k562_rese, file)
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
  theme(text = element_text(size = 16))+
  coord_cartesian(ylim = c(0, 1), xlim = c(0,1), expand = FALSE)+
  theme(legend.title=element_blank())+
  labs(y= 'True Positive Rate', x = 'False Positive Rate' )

ggsave('fig1d.png', plot = plt, width = 6, height = 4, units = 'in')

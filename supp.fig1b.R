source('auroc.R')
library(ggplot2)
library(DescTools)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/bulyk'
hepg2_rese = readPeakFile(paste(data_path, 'HepG2_all.txt', sep=''))
marks = c('H2az', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k79me2', 'H3k9ac', 'H3k9me3', 'H4k20me1')

roc_summary = data.frame()

for (mark in marks){
  file = paste(data_path, '/ucsc_hepg2/', 'wgEncodeBroadHistoneHepg2', mark, 'StdPk.broadPeak.gz', sep='')
  results = compute_roc(hepg2_rese, file)
  results$mark = mark
  roc_summary = rbind(roc_summary, results)
}


hepg2_aurocs = data.frame(mark = marks, auroc = 0)

for (mark in marks){
  mark_roc = roc_summary[which(roc_summary$mark == mark), ]
  auroc = AUC(mark_roc$FPR, mark_roc$TPR)
  hepg2_aurocs[which(hepg2_aurocs$mark == mark), 'auroc'] = auroc
}

hepg2_aurocs$names = c('H2A.Z', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1')

hepg2_aurocs$disp = paste(hepg2_aurocs$names, round(hepg2_aurocs$auroc, 3), sep = ' = ')
#hepg2_aurocs = hepg2_aurocs[rev(order(hepg2_aurocs$auroc)),]

plt = ggplot(roc_summary, aes(x = FPR, y = TPR, color = mark)) +
  scale_color_discrete(labels = hepg2_aurocs$disp)+
  geom_line(linewidth = 1.1)+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  theme(text = element_text(size = 14))+
  coord_cartesian(ylim = c(0, 1), xlim = c(0,1), expand = FALSE)+
  theme(legend.title=element_blank())+
  labs(y= 'True Positive Rate', x = 'False Positive Rate' )

ggsave('supp.fig1b.png', plot = plt, width = 6, height = 4, units = 'in')

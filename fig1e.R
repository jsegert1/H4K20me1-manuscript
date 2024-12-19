source("permutation_test.R")
library(ggplot2)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

h4K20me1_permutation = histone_permutation_test_starr(s2_silencers_granges, paste(data_path, 'H4K20me1.consensus.bed', sep = ''), numPermutations = 20000, blacklist = dm3_blacklist)

h4K20me1_permutation_DLM3 = histone_permutation_test_starr(subset(s2_silencers_granges, FIMO_DLM3 != 'none'), paste(data_path, 'H4K20me1.consensus.bed', sep = ''), numPermutations = 20000, blacklist = dm3_blacklist)

h4K20me1_permutation_SuHw = histone_permutation_test_starr(subset(s2_silencers_granges, FIMO_SuHw != 'none'), paste(data_path, 'H4K20me1.consensus.bed', sep = ''), numPermutations = 20000, blacklist = dm3_blacklist)

h4K20me1_permutation_Phaser = histone_permutation_test_starr(subset(s2_silencers_granges, FIMO_Phaser != 'none'), paste(data_path, 'H4K20me1.consensus.bed', sep = ''), numPermutations = 20000, blacklist = dm3_blacklist)

s2_H4K20me1_permutation_summary = data.frame(rbind(c('all', h4K20me1_permutation$p_val, h4K20me1_permutation$cover_ratio, h4K20me1_permutation$fold_change),
                                                   c('DLM3', h4K20me1_permutation_DLM3$p_val, h4K20me1_permutation_DLM3$cover_ratio, h4K20me1_permutation_DLM3$fold_change),
                                                   c('SuHw', h4K20me1_permutation_SuHw$p_val, h4K20me1_permutation_SuHw$cover_ratio, h4K20me1_permutation_SuHw$fold_change),
                                                   c('Phaser', h4K20me1_permutation_Phaser$p_val, h4K20me1_permutation_Phaser$cover_ratio, h4K20me1_permutation_Phaser$fold_change)))

colnames(s2_H4K20me1_permutation_summary) = c('class', 'p_val', 'cover_ratio', 'fold_change')
s2_H4K20me1_permutation_summary$p_val = as.numeric(s2_H4K20me1_permutation_summary$p_val)
s2_H4K20me1_permutation_summary$cover_ratio = as.numeric(s2_H4K20me1_permutation_summary$cover_ratio)
s2_H4K20me1_permutation_summary$fold_change = as.numeric(s2_H4K20me1_permutation_summary$fold_change)

s2_H4K20me1_permutation_summary$class = c("all\nn = 347", "DLM3\nn = 267", "SuHw\nn = 175", "Phaser\nn = 15")

plt = ggplot(s2_H4K20me1_permutation_summary, aes(x = class, y = fold_change))+
  geom_point(aes(size = cover_ratio))+
  theme_bw() + labs(x = 'Class', size = 'Cover ratio', y = 'Fold enrichment')+
  theme(text = element_text(size = 14))

ggsave('fig1e.png', plot = plt, width = 5, height = 4, units = 'in')

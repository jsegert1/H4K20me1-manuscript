library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(org.Dm.eg.db)
library(clusterProfiler)
data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

S2_H4K20me1_consensus = readPeakFile(paste(data_path, 'H4K20me1.consensus.bed', sep = ''))


S2_H4K20me1_anno = annotatePeak(S2_H4K20me1_consensus, tssRegion=c(-250, 50), TxDb=TxDb.Dmelanogaster.UCSC.dm3.ensGene, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), annoDb="org.Dm.eg.db")
h4k20me1_genes = as.data.frame(S2_H4K20me1_anno)$ENTREZID

ego = enrichGO(gene = as.vector(h4k20me1_genes), OrgDb = org.Dm.eg.db, keyType = "ENTREZID", ont = 'BP')

res = ego@result


res$qvalue_log10 = -log10(res$qvalue)
res$fg_ratio = unlist(lapply(res$GeneRatio, function(x){as.numeric(unlist(strsplit(x, '/'))[1]) / as.numeric(unlist(strsplit(x, '/'))[2])}))
res = res[rev(order(res$fg_ratio)),]

plt = ggplot(res[1:10,], aes(x = fg_ratio, y = Description, size = Count, color = qvalue_log10))+
  geom_point()+
  scale_y_discrete(limits = rev(res$Description[1:10]))+
  theme_bw()+
  labs(x = 'Fraction of H4K20me1+ genes', y = '', color = '-log10(q)')+
  scale_color_gradientn(colors = c("blue", "red"))+
  theme(text = element_text(size = 12))+
  scale_size(range = c(2, 10))+
  theme(axis.text.y = element_text(face="bold"))

ggsave("fig2f.png", plt, width = 6.5, height = 5, units = 'in')


S2_expression = read.delim(paste(paste(data_path, "Table_S1_G3-2023-404112_S2_expression.txt", sep = '')))

S2_expression$H4K20me1 = S2_expression$Gene.ID %in% data.frame(S2_H4K20me1_anno)$ENTREZID

plt2 = ggplot(S2_expression, aes(x=H4K20me1, y = S2.TPM+0.1, fill = H4K20me1))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position="none")+
  labs(x = '', y = 'S2 Transcripts per Million')+
  theme(text = element_text(size = 12))

ggsave("fig2e.png", plt2, width = 4, height = 5, units = 'in')

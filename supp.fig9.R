library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ChIPseeker)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

k562_h4k20me1_peaks = readPeakFile(paste(data_path, 'encode_chip/k562_h4k20me1.bed', sep=''))

hepg2_h4k20me1_peaks = readPeakFile(paste(data_path, 'encode_chip/hepg2_h4k20me1.bed', sep=''))

k562_H4K20me1_anno = annotatePeak(k562_h4k20me1_peaks, TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, tssRegion = c(-500,500), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic", "Downstream"))

colnames(k562_H4K20me1_anno@detailGenomicAnnotation) = c("Genic", "Intergenic", "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal intergenic")
k562_H4K20me1_anno@detailGenomicAnnotation = subset(k562_H4K20me1_anno@detailGenomicAnnotation, select = -c(Downstream))

hepg2_H4K20me1_anno = annotatePeak(hepg2_h4k20me1_peaks, TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, tssRegion = c(-500,500), genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic", "Downstream"))

colnames(hepg2_H4K20me1_anno@detailGenomicAnnotation) = c("Genic", "Intergenic", "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal intergenic")
hepg2_H4K20me1_anno@detailGenomicAnnotation = subset(hepg2_H4K20me1_anno@detailGenomicAnnotation, select = -c(Downstream))


plt = upsetplot(k562_H4K20me1_anno)
ggsave('supp.fig9a.png', plot = plt, width = 5, height = 4, units = 'in')

plt2 = upsetplot(hepg2_H4K20me1_anno)
ggsave('supp.fig9b.png', plot = plt2, width = 5, height = 4, units = 'in')


plt3 = ggplot(data.frame(k562_H4K20me1_anno), aes(x= distanceToTSS))+
  geom_histogram(bins = 50)+
  theme_bw()+
  labs(x = "Distance to TSS (bp)")+
  scale_x_continuous(limits = c(-5e5, 5e5))+
  scale_y_break(c(5000, 15000))+
  theme(
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank())+
  theme(text = element_text(size = 16))

ggsave('supp.fig9c.png', plot = plt3, width = 5, height = 4, units = 'in')

plt4 = ggplot(data.frame(hepg2_H4K20me1_anno), aes(x= distanceToTSS))+
  geom_histogram(bins = 50)+
  theme_bw()+
  labs(x = "Distance to TSS (bp)")+
  scale_x_continuous(limits = c(-5e5, 5e5))+
  scale_y_break(c(5000, 17000))+
  theme(
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank())+
  theme(text = element_text(size = 16))

ggsave('supp.fig9d.png', plot = plt4, width = 5, height = 4, units = 'in')

plt5 = plotPeakProf2(peak = k562_h4k20me1_peaks, upstream = rel(0.5), downstream = rel(0.5),
                     conf = 0.95, by = "gene", type = "body", nbin = 200,
                     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,ignore_strand = F)+theme(text = element_text(size = 16))+labs(x = 'Normalized Gene Length')

ggsave('supp.fig9e.png', plot = plt5, width = 5, height = 4, units = 'in')

plt6 = plotPeakProf2(peak = hepg2_h4k20me1_peaks, upstream = rel(0.5), downstream = rel(0.5),
                     conf = 0.95, by = "gene", type = "body", nbin = 200,
                     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,ignore_strand = F)+theme(text = element_text(size = 16))+labs(x = 'Normalized Gene Length')

ggsave('supp.fig9f.png', plot = plt6, width = 5, height = 4, units = 'in')

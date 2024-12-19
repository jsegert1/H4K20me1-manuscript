library(ggplot2)
library(ChIPseeker)
library(GenomicFeatures)
library(ggbreak)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

S2_H4K20me1_consensus = read.table(paste(data_path, 'H4K20me1.consensus.bed', sep = ''))
S2_H4K20me1_consensus$V1 = unlist(lapply(S2_H4K20me1_consensus$V1, function(x){substr(x, 4, nchar(x))})) #remove chr to match flybase
colnames(S2_H4K20me1_consensus) = c('seqnames', 'start', 'end')
S2_H4K20me1_consensus = makeGRangesFromDataFrame(S2_H4K20me1_consensus, keep.extra.columns = T)

flybase_txdb = makeTxDbFromGFF(paste(data_path, 'dmel-all-r5.57.gff.gz', sep = ''))

S2_H4K20me1_anno = annotatePeak(S2_H4K20me1_consensus, tssRegion=c(-250, 50), TxDb=flybase_txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"))

colnames(S2_H4K20me1_anno@detailGenomicAnnotation) = c("Genic", "Intergenic", "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal intergenic")
S2_H4K20me1_anno@detailGenomicAnnotation = subset(S2_H4K20me1_anno@detailGenomicAnnotation, select = -c(Downstream))

plt = upsetplot(S2_H4K20me1_anno)

ggsave('fig2c.png', plot = plt, width = 5, height = 4, units = 'in')

plt2 = ggplot(data.frame(S2_H4K20me1_anno), aes(x= distanceToTSS))+
  geom_histogram(bins = 50)+
  theme_bw()+
  labs(x = "Distance to TSS (bp)")+
  scale_x_continuous(limits = c(-1e4, 1e4))+
  scale_y_break(c(150, 2100))+
  theme(
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank())+
  theme(text = element_text(size = 16))

ggsave('supp.fig8.png', plot = plt2, width = 5, height = 4, units = 'in')


plt3 = plotPeakProf2(peak = S2_H4K20me1_consensus, upstream = rel(0.5), downstream = rel(0.5),
              conf = 0.95, by = "gene", type = "body", nbin = 200,
              TxDb = flybase_txdb,ignore_strand = F)+theme(text = element_text(size = 16))+labs(x = 'Normalized Gene Length')


ggsave('fig2d.png', plot = plt3, width = 6, height = 4, units = 'in')


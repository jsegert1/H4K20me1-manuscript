library(ggplot2)
library(ChIPseeker)
library(GenomicFeatures)
library(ggbreak)

data_path = '/data/bulyk/Julian/H4K20me1_analysis/data/'

flybase_txdb = makeTxDbFromGFF(paste(data_path, 'dmel-all-r5.57.gff.gz', sep = ''))

S2_H4K20me1_consensus = read.table(paste(data_path, 'S2_tracks/H4K20me1.consensus.bed', sep = ''))
S2_H4K20me1_consensus$V1 = unlist(lapply(S2_H4K20me1_consensus$V1, function(x){substr(x, 4, nchar(x))})) #remove chr to match flybase
colnames(S2_H4K20me1_consensus) = c('seqnames', 'start', 'end')
S2_H4K20me1_consensus = makeGRangesFromDataFrame(S2_H4K20me1_consensus, keep.extra.columns = T)


S2_H3K36me3_peaks = read.table(paste(data_path, 'Lesley_H3K36me3_dm3.bed', sep = ''))
S2_H3K36me3_peaks$V1 = unlist(lapply(S2_H3K36me3_peaks$V1, function(x){substr(x, 4, nchar(x))})) #remove chr to match flybase
colnames(S2_H3K36me3_peaks) = c('seqnames', 'start', 'end')
S2_H3K36me3_peaks = makeGRangesFromDataFrame(S2_H3K36me3_peaks, keep.extra.columns = T)


makeVennDiagram(list(S2_H4K20me1_consensus, S2_H3K36me3_peaks), 
                NameOfPeaks = c('H4K20me1', 'H3K36me3'), 
                cat.fontface = 2, 
                cat.fontfamily = 'arial', 
                sub.fontfamily = 'arial', 
                main.fontfamily = 'arial', 
                euler.d = T, 
                margin = 0.1, 
                scaled = T, 
                cat.just = list(c(.5,-1.75), c(.5,-1.75)),
                cex = 1.5,
                cat.cex = 1.5)


S2_H3K36me3_anno = annotatePeak(S2_H3K36me3_peaks, tssRegion=c(-250, 50), TxDb=flybase_txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"))

colnames(S2_H3K36me3_anno@detailGenomicAnnotation) = c("Genic", "Intergenic", "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal intergenic")
S2_H3K36me3_anno@detailGenomicAnnotation = subset(S2_H3K36me3_anno@detailGenomicAnnotation, select = -c(Downstream))

plt = ggplot(data.frame(S2_H3K36me3_anno), aes(x= distanceToTSS))+
  geom_histogram(bins = 50)+
  theme_bw()+
  labs(x = "Distance to TSS (bp)")+
  scale_x_continuous(limits = c(-1e4, 1e4))+
  scale_y_break(c(250, 1700))+
  theme(
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank())+
  theme(text = element_text(size = 16))

ggsave('supp.fig10b.png', plot = plt, width = 8, height = 4, units = 'in')

plt2 = upsetplot(S2_H4K20me1_anno)+theme(font.size = 16)
ggsave('supp.fig10c.png', plot = plt2, width = 6, height = 3, units = 'in')

plt3 = plotPeakProf(peak = list('H4K20me1' = S2_H4K20me1_consensus, 'H3K36me3' = S2_H3K36me3_peaks), upstream = rel(0.5), downstream = rel(0.5),
             conf = 0.95, by = "gene", type = "body", nbin = 200,
             TxDb = flybase_txdb,ignore_strand = F, facet = 'none')+theme(text = element_text(size = 16))+labs(x = 'Normalized Gene Length')


ggsave('supp.fig10d.png', plot = plt3, width = 10, height = 4, units = 'in')

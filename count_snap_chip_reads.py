import subprocess
import pandas as pd

fastqPath = "" #fill in your path

libs = ['Input-R1', 'H4K20me1-Abcam9051-R1', 'H4K20me1-AM-R1', 'H4K20me1-Thermo-R1', 'H3K27me3-R1', 'Input-R2', 'H4K20me1-Abcam9051-R2', 'H4K20me1-AM-R2', 'H4K20me1-Thermo-R2', 'H3K27me3-R2']

spikein_mods = ["Unmodified", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me1", "H3K9me2", "H3K9me3", "H3K27me1", "H3K27me2", "H3K27me3", "H3K36me1", "H3K36me2", "H3K36me3", "H4K20me1", "H4K20me2", "H4K20me3"]
spikein_nucs = []
for mod in spikein_mods:
    spikein_nucs.append(mod+'_A')
    spikein_nucs.append(mod+'_B')

barcodes = ["TCGCGCGTAACGACGTACCGT", "CGCGATACGACCGCGTTACGCG", "CGACGTTAACGCGTTTCGTACG", "CGCGACTATCGCGCGTAACGCG", "CCGTACGTCGTGTCGAACGACG", "CGATACGCGTTGGTACGCGTAA", "TAGTTCGCGACACCGTTCGTCG", "TCGACGCGTAAACGGTACGTCG", "TTATCGCGTCGCGACGGACGTA", "CGATCGTACGATAGCGTACCGA", "CGCATATCGCGTCGTACGACCG", "ACGTTCGACCGCGGTCGTACGA", "ACGATTCGACGATCGTCGACGA", "CGATAGTCGCGTCGCACGATCG", "CGCCGATTACGTGTCGCGCGTA", "ATCGTACCGCGCGTATCGGTCG", "CGTTCGAACGTTCGTCGACGAT", "TCGCGATTACGATGTCGCGCGA", "ACGCGAATCGTCGACGCGTATA", "CGCGATATCACTCGACGCGATA", "CGCGAAATTCGTATACGCGTCG", "CGCGATCGGTATCGGTACGCGC", "GTGATATCGCGTTAACGTCGCG", "TATCGCGCGAAACGACCGTTCG", "CCGCGCGTAATGCGCGACGTTA", "CCGCGATACGACTCGTTCGTCG", "GTCGCGAACTATCGTCGATTCG", "CCGCGCGTATAGTCCGAGCGTA", "CGATACGCCGATCGATCGTCGG", "CCGCGCGATAAGACGCGTAACG", "CGATTCGACGGTCGCGACCGTA", "TTTCGACGCGTCGATTCGGCGA"]


fastqs = []
for lib in libs:
    fastqs.append(lib+'_R1_001.fastq')
    fastqs.append(lib+'_R2_001.fastq')

spikein_counts = pd.DataFrame(index = spikein_nucs, columns = fastqs)

for i in range(len(fastqs)):
    fname = fastqPath + fastqs[i]
    for j in range(len(spikein_nucs)):
        spikein_counts.iloc[j, i] = int(subprocess.run('grep -c %s %s' % (barcodes[j], fname), stdout = subprocess.PIPE, shell=True).stdout)


spikein_counts.to_csv('spike_counts.txt', sep = '\t')
library(regioneR)
library(roxygen2)
library(GenomicRanges)

# Regulator islands (hg38)
rislands <- read.csv('/kriegsteinlab/data3/juanmoriano/MultiLayered_IndirectNeuro/data/regulatory_islands/regulatory_islands_completeINFO.tsv', header=TRUE, sep = '\t')
rislands.hg38 <- rislands[,c(5,6,7,4)] #hg3
names(rislands.hg38) <- c('chr','start','end','id')
head(rislands.hg38)
rislands.hg38 <- na.omit(rislands.hg38)
rislands.hg38 <- rislands.hg38[!duplicated(rislands.hg38),]
rislands.hg38.gr <- with(rislands.hg38, GRanges(chr, IRanges(start, end), id=id))
genome(x = rislands.hg38.gr) <- 'hg38'

# Atlas of open chromatin regions (hg38)
consensus <- read.csv('/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/consensus_signals.bed', header=FALSE, sep = '\t')
names(consensus) <- c('chr','start','end','id')
consensus.gr <- with(consensus, GRanges(chr, IRanges(start, end), id=id))
genome(x = consensus.gr) <- 'hg38'


# Are RR of genes from GEPs statistically associated to regulatory islands more often than by chance?
genes.in.ipc_late <- read.csv('/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/EarlyLate_piNMF_IP.tsv', header=TRUE, sep = '\t')
genes.in.org_late <- read.csv('/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/EarlyLate_piNMF_oRG.tsv', header=TRUE, sep = '\t')
names(genes.in.ipc_late) <- c('chr','start','end')
genes.in.ipc_late <- genes.in.ipc_late[!duplicated(genes.in.ipc_late),]
names(genes.in.org_late) <- c('chr','start','end')
genes.in.org_late <- genes.in.org_late[!duplicated(genes.in.org_late),]
genes.in.org_late.gr <- with(genes.in.org_late, GRanges(chr, IRanges(start, end)))
genes.in.ipc_late.gr <- with(genes.in.ipc_late, GRanges(chr, IRanges(start, end)))
genome(x = genes.in.org_late.gr) <- 'hg38'
genome(x = genes.in.ipc_late.gr) <- 'hg38'

#oRG
pt0 <- permTest(A=genes.in.org_late.gr, 
                ntimes=10000, 
                randomize.function=resampleRegions, 
                universe=consensus.gr,
                evaluate.function=numOverlaps, 
                B=rislands.hg38, genome='hg38',
                verbose=FALSE)
pdf(file = '/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/Permutation_RRislands_oRG.pdf', width = 5, height = 5)
plot(pt0)
dev.off()

#IP
pt1 <- permTest(A=genes.in.ipc_late.gr, 
                ntimes=10000, 
                randomize.function=resampleRegions, 
                universe=consensus.gr,
               evaluate.function=numOverlaps, 
               B=rislands.hg38, genome='hg38',
               verbose=FALSE)

pdf(file = '/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/Permutation_RRislands_IP.pdf', width = 5, height = 5)
plot(pt1)
dev.off()


# Regulator islands (hg19) - for inteserction with PosSel and Deserts of introgression
rislands.hg19 <- rislands[,c(1,2,3,4)]
names(rislands.hg19) <- c('chr','start','end','id')
head(rislands.hg19)
rislands.hg19 <- na.omit(rislands.hg19)
rislands.hg19 <- rislands.hg19[!duplicated(rislands.hg19),]
rislands.hg19.gr <- with(rislands.hg19, GRanges(chr, IRanges(start, end), id=id))
genome(x = rislands.hg19.gr) <- 'hg19'

# Regulatory islands' overlap on Pos Sel and Deserts of introgression
rpossel <- read.csv('/kriegsteinlab/data3/juanmoriano/data_CBL/BEDfiles/2020_pey_coords.bed', header=FALSE, sep = '\t')
rdeserts <- read.csv('/kriegsteinlab/data3/juanmoriano/data_CBL/BEDfiles/2020_akeydeserts_coords.bed', header=FALSE, sep = '\t')
names(rpossel) <- c('chr','start','end')
names(rdeserts) <- c('chr','start','end')
rpossel.gr <- with(rpossel, GRanges(chr, IRanges(start, end)))
rdeserts.gr <- with(rdeserts, GRanges(chr, IRanges(start, end)))
genome(x = rpossel.gr) <- 'hg19'
genome(x = rdeserts.gr) <- 'hg19'

pt2 <- permTest(A=rislands.hg19.gr, B=rpossel.gr, randomize.function=randomizeRegions,
               evaluate.function=numOverlaps, genome='hg19', ntimes=10000)

pdf(file = '/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/Permutation_RRislands_PosSel.pdf', width = 5, height = 5)
plot(pt2)
dev.off()

pt3 <- permTest(A=rislands.hg19.gr, B=rdeserts.gr, randomize.function=randomizeRegions,
               evaluate.function=numOverlaps, genome='hg19', ntimes=10000)
pdf(file = '/kriegsteinlab/data3/juanmoriano/data_CBL/indNeuro_tmp/3.Evo/Permutation_RRislands_Deserts.pdf', width = 5, height = 5)
plot(pt3)
dev.off()

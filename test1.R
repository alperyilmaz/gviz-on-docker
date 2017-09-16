##1. import data
library(Gviz)
library(rtracklayer)
library(trackViewer)

# for alper 
# setwd("/media/data/tez-seyit-chipseq/10TF_study/")

start<-48651500
end<-48658524
Ranges <- GRanges("chr6", IRanges(start,end))

data_bw <- importScore(file.path("wgEncodeCrgMapabilityAlign36mer.bigWig"),format="BigWig",ranges=Ranges)
strand(data_bw$dat) <- "-"


## rtracklayer can't read bigWig files on a Windows computer. Type in  ?`BigWigFile-class` to get the help
## ?import dedi?imizde import.bw diye bir fonksiyon var, onu kullanmama?z gerekti?ini s?yl?yor
## onu da bi denedim yine hata veriyor

data_bg <- importScore("SRR054916_MM_H0_slop_mappable.coverage.bedgraph", format="bedGraph", ranges=Ranges)
##file is too huge. Please consider to use bedtools or bedops to subset the data.hata da vermiyor ama bitmiyor da...

dat <- coverageGR(data_bg$dat)

data_bg$dat <- dat[strand(dat)=="+"]
data_bg$dat2 <- dat[strand(dat)=="-"]


##2.build gene model

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)

mmusculus <- geneModelFromTxdb(TxDb.Mmusculus.UCSC.mm9.knownGene,org.Mm.eg.db,gr=Ranges)

##3.view the tracks
viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
vp <- viewTracks(trackList(data_bw, data_bg, mmusculus), gr=Ranges, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE)
addGuideLine(c(48651584,48658224), vp=vp)
addArrowMark(list(x=48651584, y=2),label="label", col="blue",vp=vp)

#for highlighting
#idxTrack <- IdeogramTrack(genome = "mm9", chromosome = "chr6")
#axTrack <- GenomeAxisTrack()
#gtrack <- GenomeAxisTrack()
#ht <- HighlightTrack(trackList = list(data_bw, data_bg,axTrack,gtrack), start = c(48616000, 48618000), width = 7000,chromosome = 6)
#plotTracks(list(idxTrack, mmusculus, ht), from = start,
#             + to = end)


#BAM file

library(biomaRt)
bmt <- BiomartGeneRegionTrack(genome = "mm9", chromosome = "chr6", start = start, end = end, filter = list(with_ox_refseq_mrna = TRUE), stacking = "dense")
alTrack <- AlignmentsTrack(system.file(package = "Gviz","extdata", "SRR054916_MM_H0_sorted.bam"), isPaired = TRUE)
#dosya format?n? tan?mad???n? s?yl?yor..

library(BSgenome.Mmusculus.UCSC.mm9)
sTrack <- SequenceTrack(Mmusculus)

plotTracks(c(bmt, alTrack), from = start + 12700,
           + to = start + 15200, chromosome = "chr6", reverseStacking = TRUE,
           + col.mates = "purple", col.gap = "orange", type = "pileup")

plotTracks(c(alTrack, sTrack), chromosome = "chr6",
           + from = start, to = end, cex = 0.5, min.height = 8)

sTrack <- SequenceTrack(mmusculus)
plotTracks(mmusculus)

data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = "hg19",
                           chromosome = "chr7", name = "foo")
plotTracks(grtrack)

ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
plotTracks(ideoTrack, from = 8.5e+07, to = 1.29e+08)

#### the part below is working

# from https://support.bioconductor.org/p/58768/
library(Gviz)

chr <- "chr12"
start<-80695688 		
end<-80695877
Ranges <- GRanges(chr, IRanges(start,end))

bgFile <- "SRR054916_MM_H0_slop_mappable.coverage.bedgraph"
dTrack2 <- DataTrack(range = bgFile, genome = "mm9",
                     type = "mountain", chromosome = chr, name = "bedGraph", ylim=c(0,500))

gtrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "mm9", chromosome = chr)
library(biomaRt)
bmt <- BiomartGeneRegionTrack(genome = "mm9", chromosome = chr, start = start, end = end, filter = list(with_ox_refseq_mrna = TRUE), stacking = "dense")
alTrack <- AlignmentsTrack("SRR054916_MM_H0_sorted.bam", chromosome = chr, start = start, end = end, isPaired = TRUE)

mappability <- DataTrack("wgEncodeCrgMapabilityAlign36mer.bigWig", strand="+",
                         genome="mm9", name="mappability", type="histogram")

conservation <- UcscTrack(genome = "mm9", chromosome = chr,
                          track = "Conservation", table = "phyloP30wayPlacental",
                          from = start, to = end, trackType = "DataTrack",
                          start = "start", end = "end", data = "score",
                          type = "hist", window = "auto", col.histogram = "darkblue",
                          fill.histogram = "darkblue", ylim = c(-3.7, 4),
                          name = "Conservation")

gcContent <- UcscTrack(genome = "mm9", chromosome = chr,
                       track = "GC Percent", table = "gc5Base", from = start,
                       to = end, trackType = "DataTrack", start = "start",
                       end = "end", data = "score", type = "hist", window = -1,
                       windowSize = 1500, fill.histogram = "black",
                       col.histogram = "black", ylim = c(30, 70), name = "GC Percent")
pdf('chr12-80mb.pdf',width=15,height=6)
plotTracks(list(ideoTrack, gtrack, bmt, dTrack2, alTrack, mappability, conservation, gcContent), from = start-500, to = end+500,
            chromosome=chr)
dev.off() #Clost pdf after plotting


## trying to get gene models from a file (or library) so that we dont depend on internet connection 
## to retrieve genes from biomart
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
genetrack <- GeneRegionTrack(txdb, genome = "mm9",
                chromosome = "chr6", name = "UCSC known genes", symbol = TRUE, showId = TRUE,background.title ="#042E8A",col="#042E8A",fill="#042E8A",fontsize=16)
plotTracks(list(ideoTrack,gtrack,genetrack),from=20000000,to=24000000,chromosome = "chr6") #,transcriptAnnotation = "symbol")



#Create a GViz track showing how the annotated transcripts sit in the genome
gtrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "mm9", chromosome = "chr16")
#options(ucscChromosomeNames=FALSE)
genetrack <- GeneRegionTrack("Mus_musculus.GRCm38.90.modified.gtf.gz",chromsome="chr16",name="Transcripts",transcriptAnnotation="transcript", fontsize=16,thinBoxFeature=c("five_prime_utr,three_prime_utr"),collapse=TRUE) #background.title ="#619CFF",fill="#619CFF",col="black",size=0.5,
genetrack <- BiomartGeneRegionTrack(genome = "mm9", chromosome = "chr16", start = 56686000, end = 56693000, filter = list(with_ox_refseq_mrna = TRUE),transcriptAnnotation="transcript",background.title ="#619CFF",fill="#619CFF",col="black",size=0.5)

#Open the R pdf saving functions
pdf('coverage_ref2.pdf',width=15,height=6)
plotTracks(list(ideoTrack,gtrack, genetrack),from=56686000,to=56693000,chromosome="chr16") #,type=c("coverage"))
dev.off() #Clost pdf after plotting

## As a result, gtf solution cannot render UTR correctly.. check coverage_ref1.pdf and coverage_ref2.pdf files for comparison.
## thus biomart is the only correct solution

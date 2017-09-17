#### the part below is working

# from https://support.bioconductor.org/p/58768/
library(Gviz)
library(biomaRt)

chr <- "chr6"
start<-48617000
end<-48620000
Ranges <- GRanges(chr, IRanges(start,end))

bgFile <- "regions-for-publication-coverage.bedgraph"
dTrack2 <- DataTrack(range = bgFile, 
                     genome = "mm9",
                     type = "mountain", 
                     chromosome = chr, 
                     name = "bedGraph", 
                     ylim=c(0,100))

gtrack <- GenomeAxisTrack()

ideoTrack <- IdeogramTrack(genome = "mm9", chromosome = chr)

bmt <- BiomartGeneRegionTrack(genome = "mm9", 
                              chromosome = chr, 
                              start = start, 
                              end = end, 
                              filter = list(with_ox_refseq_mrna = TRUE), 
                              stacking = "dense")

alTrack <- AlignmentsTrack("regions-for-publication-sorted.bam", 
                           chromosome = chr, 
                           start = start, 
                           end = end, 
                           isPaired = TRUE)

mappability <- DataTrack("regions-for-publication-mappability.bw",
                         genome="mm9", 
                         name="mappability", 
                         type="histogram")  # strand="+"

conservation <- UcscTrack(genome = "mm9", chromosome = chr,
                          track = "Conservation", table = "phyloP30wayPlacental",
                          from = start, to = end, trackType = "DataTrack",
                          start = "start", end = "end", data = "score",
                          type = "hist", window = "auto", col.histogram = "darkblue",
                          fill.histogram = "darkblue", ylim = c(-3.7, 4),
                          name = "Conservation")

# gcContent <- UcscTrack(genome = "mm9", chromosome = chr,
#                        track = "GC Percent", table = "gc5Base", from = start,
#                        to = end, trackType = "DataTrack", start = "start",
#                        end = "end", data = "score", type = "hist", window = -1,
#                        windowSize = 1500, fill.histogram = "black",
#                        col.histogram = "black", ylim = c(30, 70), name = "GC Percent")

#pdf('chr12-80mb.pdf',width=15,height=6)
plotTracks(list(ideoTrack, gtrack, bmt, dTrack2, alTrack, mappability, conservation), from = start+450, to = end-700,
            chromosome=chr, background.title ="#042E8A",fontsize=14)
#dev.off() #Clost pdf after plotting


## trying to get gene models from a file (or library) so that we dont depend on internet connection 
## to retrieve genes from biomart. I tried TxDb.Mmusculus.UCSC.mm9.knownGene or downloading gtf file and using it locally.
## As a result, gtf solution cannot render UTR correctly.. check coverage_ref1.pdf and coverage_ref2.pdf files for comparison.
## thus biomart is the only viable solution
## TxDb should be working too
# library(TxDb.Mmusculus.UCSC.mm9.knownGene)
# gtTrack <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm9.knownGene, chromosome="chr16", start=16858904, end=16895526)
# plotTracks(gtTrack,transcriptAnnotation="name")

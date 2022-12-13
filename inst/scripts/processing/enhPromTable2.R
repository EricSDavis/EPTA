## Load required libraries
library(data.table)
library(GenomicRanges)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(hictoolsr)

## Shorten txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Input two sets of GRanges for creating interaction pairs ------------------------------

## Define promoter regions from differential rna-seq data
prom <-
  "inst/extdata/rna/LIMA_RNA_diffGenes_binned.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb)) |>
  promoters(upstream = 2000, downstream = 200)

colnames(mcols(prom)) <- paste0("RNA_", colnames(mcols(prom)))


## Define enhancers as regulatory regions - promoters
enh <-
  "inst/extdata/atac_h3k27ac/LIMA_H3K27ac_diffPeaks_overlap.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb)) |>
  subsetByOverlaps(ranges = prom,
                   invert = TRUE)

colnames(mcols(enh)) <- paste0("K27_", colnames(mcols(enh)))


## Calculate interaction pairs between enh and prom within 2Mb
ov <- findOverlaps(enh, prom, maxgap = 2e6)
ep <- GInteractions(anchor1 = enh[queryHits(ov)],
                    anchor2 = prom[subjectHits(ov)])


## Annotate interaction pairs with additional features -----------------------------------

## Read in loops
loops <-
  "inst/extdata/hic/LIMA_HiC_diffLoops.tsv" |>
  fread() |>
  as_ginteractions(starts.in.df.are.0based = TRUE,
                   keep.extra.columns = TRUE)

colnames(mcols(loops)) <- paste0("HIC_", colnames(mcols(loops)))


## Add loop data (# slow :( ~ 12 seconds)
loopOV <- findOverlaps(query = ep, subject = loops)
system.time({
  mcols(ep[queryHits(loopOV),]) <-
    cbind(mcols(ep[queryHits(loopOV),]),
          mcols(loops[subjectHits(loopOV),]))
})


## Put EP pairs into hic bins
binnedEP <- binBedpe(bedpe = ep, res = 10e3, a1Pos = 'center', a2Pos = 2000)

## Throw out interactions that bin past the end of a chromosome
binnedEP <-
  lapply(c('first', 'second'),
         \(x) which(width(anchors(binnedEP, x)) != 10e3+1)) |>
  unlist() |>
  unique() |>
  {\(x) binnedEP[-x]}()

system.time({
  epCounts <-
    BiocParallel::bplapply(X = c(1:22, 'X', 'Y'), FUN = \(x) {
      message("chr",x)
      extractCounts(bedpe = binnedEP,
                    hic = list.files("../../LIMA/data/storage/hic/hic_maps/",
                                     full.names = TRUE),
                    chroms = x,
                    res = 10e3,
                    norm = 'KR',
                    matrix = 'observed')
    })

  epCounts <- do.call(c, epCounts)
})


save(epCounts, file = "data/epCounts.rda")

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
genes <-
  "inst/extdata/rna/LIMA_RNA_diffGenes_binned.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb))

colnames(mcols(genes)) <- paste0("RNA_", colnames(mcols(genes)))

## Define enhancers as regulatory regions - promoters
enh <-
  "inst/extdata/atac_h3k27ac/LIMA_H3K27ac_diffPeaks_overlap.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb))

colnames(mcols(enh)) <- paste0("K27_", colnames(mcols(enh)))

## Compute overlaps between promoters and enhancers
ov <- findOverlaps(query = promoters(x = genes, upstream = 2000, downstream = 200),
                   subject = enh)

## Calculate interaction pairs between enh and prom within 1Mb
pgPairs <- GInteractions(anchor1 = enh[subjectHits(ov)],
                         anchor2 = genes[queryHits(ov)])

## Save object
save(pgPairs, file = "data/pgPairs.rda")

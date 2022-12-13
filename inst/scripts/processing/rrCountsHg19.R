## Load required libraries
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Read in regulatory regions with atac and k27 counts
k27 <- system.file("extdata/h3k27ac/counts/atacPeaks_H3K27AcCounts_500bp.txt",
                   package = 'EPTA') |> fread()
atac <- system.file("extdata/atac/counts/atacPeaks_atacCounts.txt",
                    package = 'EPTA') |> fread()


## Combine atac and k27 counts & convert to GRanges
rr <- makeGRangesFromDataFrame(df = merge(atac, k27),
                               keep.extra.columns = TRUE,
                               seqinfo = seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))

## Rename and save rna table
rrCountsHg19 <- rr
save(rrCountsHg19, file = "data/rrCountsHg19.rda")

## Load packages
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(tidyr)
library(purrr)
library(plotgardener)

## Load gene and enh locations -----------------------------------------------------------

## Shorten txdb name
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Differential rna-seq data
genes <-
  "inst/extdata/rna/LIMA_RNA_diffGenes_binned.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb)) |>
  resize(width = 1)
colnames(mcols(genes)) <- paste0("RNA_", colnames(mcols(genes)))

## Define promoters
prom <- promoters(x = genes, upstream = 2000, downstream = 200)

## Define enhancers as regulatory regions - promoters
enh <-
  "inst/extdata/atac_h3k27ac/LIMA_H3K27ac_diffPeaks_overlap.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb)) |>
  subsetByOverlaps(ranges = prom,
                   invert = TRUE)
colnames(mcols(enh)) <- paste0("K27_", colnames(mcols(enh)))


## Load filtered epCounts, abc, and d metrics --------------------------------------------
load("data/abcdResults.rda", verbose = TRUE)

## Add ABCD to epCounts
epABCD <- epCounts
mcols(epABCD) <- cbind(mcols(epABCD), abc)
epABCD$maxABC <- rowMaxs(as.matrix(mcols(epABCD)[grep("ABC", colnames(mcols(epABCD)))]))
epABCD$dynamics <- 1-d

## Correct enhancer and gene locations in epCounts object
gi <- GInteractions(anchor1 = granges(enh[match(epABCD$anchor1.K27_name,
                                                enh$K27_name)]),
                    anchor2 = granges(genes[match(epABCD$anchor2.RNA_ensembl,
                                                  genes$RNA_ensembl)]))
mcols(gi) <- mcols(epABCD)
epABCD <- gi


## Visualize maxABC vs Dynamics for a single gene ----------------------------------------

## Define min/max normalization function
sc <- \(x){rm <- rowMins(x); (x-rm)/(rowMaxs(x)-rm)} #matrix-wide



## Extract and normalize abc and rna
abc <- sc(as.matrix(goi[,grep("^ABC", colnames(goi))]))
rna <- sc(as.matrix(goi[,grep("RNA", colnames(goi))]))

## Plot points of interest
plotScores <- \(x) {
  tibble(abc = abc[x,],
         rna = rna[x,],
         timepoint = 1:8) |>
    pivot_longer(cols = c(abc, rna), names_to = "type") |>
    ggplot(aes(x = timepoint, y = value, col = type)) +
    scale_color_manual(values = c("black", "red")) +
    scale_x_continuous(breaks = 1:8) +
    geom_line() +
    geom_point(fill = "white", shape = 21, size = 2, stroke = 1) +
    labs(y = "",
         x = "Timepoint")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank())
}





## Visualize this region -----------------------------------------------------------------

## Define parameters
p <- pgParams(gene = "IL1B",
              geneBuffer = c(80e3, 120e3),
              assembly = 'hg19',
              linecolor = NA,
              width = 9,
              length = 2,
              height = 2,
              x = 0.5,
              resolution = 5e3,
              zrange = c(0, 1000),
              norm = 'KR')

## Select pairs to IL1B
pairs <- epABCD[epABCD$anchor2.RNA_hgnc == "IL1B",
                grep(pattern = "K27_cluster|RNA.*VST|ABC|dynamics",
                     x = colnames(mcols(epABCD)))]

## Filter pairs within region
region <- GRanges(seqnames = p$chrom,
                  ranges = IRanges(start = p$chromstart,
                                   end = p$chromend),
                  seqinfo = seqinfo(pairs))
pairs <- pairs[countOverlaps(anchors(pairs, 'first'), region, type = 'within') > 0]

# goi <- as.data.frame(mcols(pairs))

## Load ATAC peaks
atacPeaks <-
  fread("../../LIMA/data/LIMA_share/atac/peaks/LIMA_ATAC_THP1_WT_LPIF_S_MACS2_peakMerge.bed") |>
  `colnames<-`(value = c("seqnames", "start", "end")) |>
  makeGRangesFromDataFrame()


## Plot score scatterplots ---------------------------------------------------------------

center = mid(anchors(pairs)$second)[1]
plotStart = center - p$geneBuffer[1]
plotEnd = center + p$geneBuffer[2]
p1 = 113596838
p2 = 113637304
p3 = 113682586 #113222929
df <- data.frame(start = start(anchors(pairs, "first")),
                 abcScore = pairs$maxABC,
                 dynamics = log2(pairs$dynamics/0.7),
                 abcdScore  = pairs$maxABC * log2(pairs$dynamics/0.7))
par(mfrow = c(3,1), pch = 19)
plot(df$start, df$abcScore, main = "max ABC", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(p1, df$abcScore[df$start == p1], col = "blue")
points(p2, df$abcScore[df$start == p2], col = "red")
points(p3, df$abcScore[df$start == p3], col = "purple")
#abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
plot(df$start, df$dynamics, main = "Dynamics", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(p1, df$dynamics[df$start == p1], col = "blue")
points(p2, df$dynamics[df$start == p2], col = "red")
points(p3, df$dynamics[df$start == p3], col = "purple")
# abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
abline(h = 0, lty = 2)
plot(df$start, df$abcdScore, main = "ABCD", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(p1, df$abcdScore[df$start == p1], col = "blue")
points(p2, df$abcdScore[df$start == p2], col = "red")
points(p3, df$abcdScore[df$start == p3], col = "purple")
# abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
abline(h = 0, lty = 2)
par(mfrow = c(1,1), pch = 1)



## Version 1 -----------------------------------------------------------------------------

## Create page
pageCreate(width = 16, height = 9)

## Plot Hic
hic <-
  plotHicRectangle(data = "/Volumes/KatieEHD/Data/LIMA_share/hic/hic0_maps/LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic",
                   params = p,
                   height = 3,
                   y = 0.5)

## Plot Signal tracks
plotSignal(data = "../../LIMA/data/LIMA_share/h3k27ac/signal/MERGE_LIMA_ChIP_h3k27ac_THP1_WT_LPIF_0120_S.bw",
           params = p,
           y=3.6,
           height = 0.9,
           range = c(0, 1000),
           fill = "firebrick")

## Plot ATAC peaks
plotRanges(data = atacPeaks,
           params = p,
           y = 4.5,
           height = 0.1,
           collapse = TRUE,
           fill = "grey60")

## Gene track
plotGenes(params = p,
          y = 4.6,
          height = 1,
          geneHighlights = data.frame(gene="IL1B", color = "orange"),
          geneBackground = "grey")

## Version 2 -----------------------------------------------------------------------------

## Create page
pageCreate(width = 16, height = 9)

## Plot Hic
hic <-
  plotHicRectangle(data = "/Volumes/KatieEHD/Data/LIMA_share/hic/hic0_maps/LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic",
                   params = p,
                   height = 3,
                   y = 0.5)

## Plot Signal tracks
plotSignal(data = "../../LIMA/data/LIMA_share/h3k27ac/signal/MERGE_LIMA_ChIP_h3k27ac_THP1_WT_LPIF_0120_S.bw",
           params = p,
           y=3.6,
           height = 0.9,
           range = c(0, 1000),
           fill = "firebrick")

## Plot ATAC peaks
plotRanges(data = atacPeaks,
           params = p,
           y = 4.5,
           height = 0.1,
           collapse = TRUE,
           linecolor = "grey60",
           fill = "grey60")

## Draw pair arches
plotPairsArches(data = swapAnchors(pairs),
                params = p,
                height = 1,
                y = 4.6,
                flip = TRUE,
                archHeight = 1,
                fill = alpha("grey60", 0.9),
                linecolor = alpha("grey60", 0.25),
                alpha = 1)

## Highlight pair of interest
plotPairsArches(data = swapAnchors(pairs)[df$start %in% c(p1, p2, p3)],
                params = p,
                height = 1,
                y = 4.6,
                flip = TRUE,
                archHeight = 1,
                fill = "blue",
                linecolor = NA,
                alpha = 1)


library(grid)

plotXY <- \(x, y, width, height, default.units = "inches", just = c("left", "top")) {
  pushViewport(viewport(x = x, y = y, width = width, height = height,
                        default.units = default.units, just = just,
                        xscale = c(p$chromstart, p$chromend), yscale = range(df$abcScore)))
  grid.rect()
  grid.xaxis()
  grid.yaxis()
  grid.points(x = df$start, y = df$abcScore, pch = 19)
  upViewport(n=1)
}

plotXY(x = 0.5, y = 9-6, width = p$width, height = 2)

## Notes ---------------------------------------------------------------------------------
## Genomics
## ABC would predict this variant
## Look at time correlation (not great for ABC selected one)
## Dynamics has some spurious correlations though
## Combining into ABCD selects the "best" candidates
## Show genomics of the connection
## This improves assignments by ABC - model for all loci

## Zoom in the size of the plot a bit

## Highlight point of interest
# plotgardener::annoHighlight(plot = hic, chrom = hic$chrom,
#                             chromstart = 113637304-1e03, chromend = 113637304+1e03, y = 0.5,
#                             height = 8)

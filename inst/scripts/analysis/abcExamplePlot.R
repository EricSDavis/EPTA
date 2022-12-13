## Load packages
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(tidyr)
library(purrr)

## Load gene and enh locations -----------------------------------------------------------

## Shorten txdb name
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Differential rna-seq data
genes <-
  "inst/extdata/rna/LIMA_RNA_diffGenes_binned.tsv" |>
  fread() |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb))
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

## Find a differential gene
colnames(mcols(epABCD))
table(epABCD$anchor2.RNA_cluster)

## Find the genes with the most joined enhancers
sort(table(epABCD$anchor2.RNA_hgnc), decreasing = TRUE) |> head(n=40)

## Visualize maxABC vs Dynamics for a single gene ----------------------------------------

## Define min/max normalization function
sc <- \(x){rm <- rowMins(x); (x-rm)/(rowMaxs(x)-rm)} #matrix-wide

## Select gene of interest
pairs <- epABCD[epABCD$anchor2.RNA_hgnc == "CITED4", #CITED4
                grep("K27_cluster|RNA.*VST|ABC|dynamics", colnames(mcols(epABCD)))]
goi <- as.data.frame(mcols(pairs))


## Visualize all putative enhancers for goi
ggplot(data = goi,
       mapping = aes(x = maxABC, y = dynamics)) +
  geom_point() +
  labs(title = 'Putative Regulators of "CITED4"',
       x = "Max ABC Score Per Timepoint",
       y = "Time Dynamics") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
# ggsave("inst/plots/abcExamplePlot_CITED4.pdf", width = 6, height = 5)

## Identify pairs of interest from each quadrant
## CITED4
urp <- which(goi$maxABC > 40 & goi$dynamics > 0.8)
lrp <- which(goi$maxABC > 35 & goi$dynamics < 0.43)
llp <- which(goi$maxABC < 15 & goi$dynamics < 0.43)
# # FOXO6
# urp <- which(goi$maxABC > 32 & goi$dynamics > 0.85)
# lrp <- which.max(goi$maxABC)
# llp <- which.min(goi$dynamics)


## Color points of interest
plot1 <-
  ggplot(data = goi,
       mapping = aes(x = maxABC, y = dynamics)) +
  geom_point(color = 'grey') +
  labs(title = 'Putative Regulators of "CITED4"',
       x = "Max ABC Score Per Timepoint",
       y = "Time Dynamics") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
plot1
# ggsave("inst/plots/abcExamplePlot_CITED4_grey.pdf", width = 6, height = 5)

plot1 +
  geom_point(data = goi[llp,], color = 'maroon', size = 3)
# ggsave("inst/plots/abcExamplePlot_CITED4_grey_llp.pdf", width = 6, height = 5)

plot1 +
  geom_point(data = goi[llp,], color = 'maroon', size = 3) +
  geom_point(data = goi[lrp,], color = 'blue', size = 3)
# ggsave("inst/plots/abcExamplePlot_CITED4_grey_lrp.pdf", width = 6, height = 5)

plot1 +
  geom_point(data = goi[llp,], color = 'maroon', size = 3) +
  geom_point(data = goi[lrp,], color = 'blue', size = 3) +
  geom_point(data = goi[urp,], color = 'purple', size = 3)
# ggsave("inst/plots/abcExamplePlot_CITED4_grey_urp.pdf", width = 6, height = 5)


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

plotScores(x = urp)
# ggsave("inst/plots/abcExamplePlot_CITED4_lines_urp.pdf", width = 6, height = 5)

plotScores(x = lrp)
# ggsave("inst/plots/abcExamplePlot_CITED4_lines_lrp.pdf", width = 6, height = 5)

plotScores(x = llp)
# ggsave("inst/plots/abcExamplePlot_CITED4_lines_llp.pdf", width = 6, height = 5)


## Visualize this region -----------------------------------------------------------------

library(plotgardener)


## Define region in parameters
p <- pgParams(gene = "CITED4",
              geneBuffer = c(2.5e3, 20e3),#1.1e06,
              assembly = 'hg19',
              width = 4,
              length = 4,
              height = 0.5,
              x = 0.5)

chipParams <- pgParams(fill = "#0063B2FF",
                       linecolor = "#0063B2FF",
                       assembly = "hg19")
rnaParams <- pgParams(fill = "#9CC3D5FF",
                      linecolor = "#9CC3D5FF",
                      assembly = "hg19")


## Init page
pageCreate(width = 5, height = 6.5, showGuides = FALSE, xgrid = 0, ygrid = 0)

## Signal track files for all timepoints
chipFiles <- list.files(path = "../../LIMA/data/LIMA_share/h3k27ac/signal",
                        pattern = ".bw",
                        full.names = TRUE)
rnaFiles <- list.files(path = "../../LIMA/data/LIMA_share/rna/signal",
                       pattern = ".bw",
                       full.names = TRUE)

## Read in data and find max yrange (use custom range)
chip <- lapply(chipFiles, \(x) readBigwig(file = x, params = p))
rna <- lapply(rnaFiles, \(x) readBigwig(file = x, params = p))

## Define ypos
ypos <- c(1, sapply(1:7, \(x) 1+(0.5+0.1)*x))

## Visualize
plotGenes(params = p, y = 0.25, height = 0.75)

pmap(list(chip, ypos), \(x,y) {
  plotSignal(data = x, params = c(p, chipParams), y = y, range = c(0,3000))
})
pmap(list(rna, ypos), \(x,y) {
  plotSignal(data = x, params = c(p, rnaParams), y = y, range = c(0,15000))
})

plotPairsArches(data = pairs, params = p,
                y = ypos[8]+p$height, flip = TRUE,
                fill = "lightgrey", linecolor = "lightgrey")
plotPairsArches(data = pairs[urp,], params = p,
                y = ypos[8]+p$height, flip = TRUE,
                fill = chipParams$fill, linecolor = chipParams$fill)

plotGenomeLabel(params = p, y = ypos[8]+p$height, scale = 'Mb')


#40320000-42330000
#pairs[!pairs$anchor1.K27_cluster %in% c('static', 'none')]


## Find the genes with the most joined enhancers
sort(table(epABCD$anchor2.RNA_hgnc), decreasing = FALSE) |> head(n=40)

pairs <- epABCD[epABCD$anchor2.RNA_hgnc == "CHRFAM7A",
                grep("K27_cluster|RNA.*VST|ABC|dynamics", colnames(mcols(epABCD)))]
goi <- as.data.frame(mcols(pairs))

## Plot whole region

pageCreate(width = 20, height = 8)
p <- pgParams(gene = "CITED4",
              geneBuffer = 50e3,#1.1e06,
              assembly = 'hg19',
              linecolor = NA,
              width = 19,
              length = 7,
              height = 4,
              x = 0.5)

chipFiles <- list.files(path = "../../LIMA/data/LIMA_share/h3k27ac/signal",
                        pattern = ".bw",
                        full.names = TRUE)
atacFiles <- list.files(path = "../../LIMA/data/LIMA_share/atac/signal",
                        pattern = ".bw",
                        full.names = TRUE)

chip <- lapply(chipFiles, \(x) readBigwig(file = x, params = p))
atac <- lapply(atacFiles, \(x) readBigwig(file = x, params = p))

plotSignal(data = chip[[1]], params = p, y=0.5, range = c(0, 1000),
           fill = alpha("#37a7db", 0.5))
plotSignal(data = atac[[1]], params = p, y=0.5, range = c(0, 1000),
           fill = alpha("indianred2", 0.5))

plotPairsArches(data = pairs[pairs$anchor1.K27_cluster == 'none',], params = p, height = 2,
                y = 4.5, flip = TRUE,
                fill = "lightgrey", linecolor = "lightgrey")
plotPairsArches(data = pairs[pairs$anchor1.K27_cluster != 'none',], params = p, height = 2,
                y = 4.5, flip = TRUE)


## With hic
pairs <- epABCD[epABCD$anchor2.RNA_hgnc == "IL1B",
                grep("K27_cluster|RNA.*VST|ABC|dynamics", colnames(mcols(epABCD)))]
goi <- as.data.frame(mcols(pairs))

pageCreate(width = 20, height = 11)
p <- pgParams(gene = "IL1B",
              geneBuffer = 200e3,#1.1e06,
              assembly = 'hg19',
              linecolor = NA,
              width = 19,
              length = 2,
              height = 2,
              x = 0.5)

chipFiles <- list.files(path = "../../LIMA/data/LIMA_share/h3k27ac/signal",
                        pattern = ".bw",
                        full.names = TRUE)
atacFiles <- list.files(path = "../../LIMA/data/LIMA_share/atac/signal",
                        pattern = ".bw",
                        full.names = TRUE)

chip <- lapply(chipFiles, \(x) readBigwig(file = x, params = p))
atac <- lapply(atacFiles, \(x) readBigwig(file = x, params = p))

hic <-
plotHicRectangle(data = "../../LIMA/data/storage/hic/hic_maps/LIMA_THP1_WT_LPIF_omega_S_0.0.0_megaMap_inter_30.hic",
                 resolution = 5e3, zrange = c(0, 1000), norm = 'KR', params = p,
                 height = 5, y = 0.5)

plotSignal(data = chip[[1]], params = p, y=5.5, range = c(0, 1000),
           fill = alpha("#37a7db", 0.5))
plotSignal(data = atac[[1]], params = p, y=5.5, range = c(0, 1000),
           fill = alpha("indianred2", 0.5))

plotgardener::annoHighlight(plot = hic, chrom = hic$chrom,
                            chromstart = 113637304-1e03, chromend = 113637304+1e03, y = 0.5,
                            height = 8)

## Make "fake" manhattan plots for abc, d, and abcd
mm <- function(x){(x-min(x))/(max(x)-min(x))}
df_abc <- data.frame(
  chrom = as.character(seqnames(anchors(pairs, 'first'))),
  pos = start(anchors(pairs, 'first')),
  p = 10^-mm(pairs$ABC_0000))

df_d <- df_abc
df_d$p <- 10^-mm(pairs$dynamics/0.7)

df_abcd <- df_abc
df_abcd$p <- 10^-mm(pairs$ABC_0000 * pairs$dynamics)


# mh <- plotManhattan(data = df_abc, params = p, y = 7.6, height = 1, sigCol = "red", sigVal = 1, range = c(0, 1.1))
# plotManhattan(data = df_d, params = p, y = 7.6, height = 1, sigCol = "orange2", sigVal = 1, range = c(0, 1.1))
# plotManhattan(data = df_abcd, params = p, y = 7.6, height = 1, sigCol = "blue", sigVal = 1, range = c(0, 1.1))
#
# annoXaxis(plot = mh, axisLine = TRUE)
# annoYaxis(plot = mh, axisLine = TRUE)

mh1 <- plotManhattan(data = df_abc, params = p, y = 7.6, height = 1, sigCol = "red", sigVal = 1, range = c(0, 1.1))
mh2 <- plotManhattan(data = df_d, params = p, y = 8.7, height = 1, sigCol = "orange2", sigVal = 1, range = c(0, 1.1))
mh3 <- plotManhattan(data = df_abcd, params = p, y = 9.8, height = 1, sigCol = "blue", sigVal = 1, range = c(0, 1.1))

annoXaxis(plot = mh1, axisLine = TRUE)
annoYaxis(plot = mh1, axisLine = TRUE)

annoXaxis(plot = mh2, axisLine = TRUE)
annoYaxis(plot = mh2, axisLine = TRUE)

annoXaxis(plot = mh3, axisLine = TRUE)
annoYaxis(plot = mh3, axisLine = TRUE)



# plotPairsArches(data = pairs, params = p, height = 2,
#                 y = 5.5, flip = TRUE, archHeight = pairs$dynamics)
# plotPairsArches(data = swapAnchors(pairs), params = p, height = 2,
#                 y = 7.5, flip = TRUE, archHeight = swapAnchors(pairs)$ABC_0000 - min(swapAnchors(pairs)$ABC_0000))

# plotPairsArches(data = swapAnchors(pairs)[pairdist(swapAnchors(pairs)) <= p$geneBuffer,], params = p, height = 2,
#                 y = 7.5, flip = TRUE, archHeight = (swapAnchors(pairs)$ABC_0000 - min(swapAnchors(pairs)$ABC_0000))[pairdist(swapAnchors(pairs)) <= p$geneBuffer[1]])
#
# (swapAnchors(pairs)$dynamics - min(swapAnchors(pairs)$dynamics))[pairdist(swapAnchors(pairs)) <= p$geneBuffer[1]]*
# (swapAnchors(pairs)$ABC_0000 - min(swapAnchors(pairs)$ABC_0000))[pairdist(swapAnchors(pairs)) <= p$geneBuffer[1]]


## Cant clip
## need to swap anchors

## 100 genes that change the most
## Arch plots

library(ggplot2)
data.frame(distance = pairdist(epABCD)[pairdist(epABCD) <= 1e06],
     dynamics = epABCD$dynamics[pairdist(epABCD) <= 1e06]) |>
  ggplot(aes(x = distance, y = dynamics)) +
  # geom_point() +
  geom_hex()
  # geom_smooth()

plot(density(epABCD$dynamics[pairdist(epABCD) <= 1e06]))
abline()

epABCD$dynamics[pairdist(epABCD) <= 1e06]

data.frame(distance = pairdist(pairs),
           dynamics = pairs$dynamics) |>
lm(formula = distance ~ dynamics)
plot(pairdist(pairs), pairs$dynamics)

lines(spline(x,y,n = 10))

epABCD


data.frame(start = start(anchors(pairs, "first")),
           dynamics = log2(pairs$dynamics/0.7)) |>
  ggplot(aes(x = start, y = dynamics)) +
  geom_point()+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 113587328)

data.frame(start = start(anchors(pairs, "first")),
           abcScore = pairs$ABC_0000) |>
  ggplot(aes(x = start, y = abcScore)) +
  geom_point()+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 113587328)

data.frame(start = start(anchors(pairs, "first")),
           abcScore = pairs$ABC_0000,
           dynamics = log2(pairs$dynamics/0.7),
           abcdScore  = pairs$ABC_0000 * log2(pairs$dynamics/0.7)) |>
  tidyr::pivot_longer(cols = c(abcScore, dynamics, abcdScore), names_to = "score") |>
  ggplot(aes(x = start, y = value)) +
  facet_grid() +
  geom_point()+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 113587328)


pairs <- epABCD[epABCD$anchor2.RNA_hgnc == "IL1B",
                grep("K27_cluster|RNA.*VST|ABC|dynamics", colnames(mcols(epABCD)))]
goi <- as.data.frame(mcols(pairs))

pairs

center = mid(anchors(pairs)$second)[1]
plotStart = center - 500000
plotEnd = center + 500000
p1 = 113596838
p2 = 113637304
p3 = 113222929
df <- data.frame(start = start(anchors(pairs, "first")),
           abcScore = pairs$maxABC,
           dynamics = log2(pairs$dynamics/0.7),
           abcdScore  = pairs$maxABC * log2(pairs$dynamics/0.7))
par(mfrow = c(3,1), pch = 19)
plot(df$start, df$abcScore, main = "max ABC", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(df$start[which.max(df$abcdScore)], df$abcScore[which.max(df$abcdScore)], col = "blue")
points(df$start[which.max(df$abcScore)], df$abcScore[which.max(df$abcScore)], col = "red")
points(df$start[38], df$abcScore[38], col = "purple") #which.max(df$dynamics)
#abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
plot(df$start, df$dynamics, main = "Dynamics", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(df$start[which.max(df$abcdScore)], df$dynamics[which.max(df$abcdScore)], col = "blue")
points(df$start[which.max(df$abcScore)], df$dynamics[which.max(df$abcScore)], col = "red")
points(df$start[38], df$dynamics[38], col = "purple")
# abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
abline(h = 0, lty = 2)
plot(df$start, df$abcdScore, main = "ABCD", xaxs = 'i', xlim=c(plotStart, plotEnd))
points(df$start[which.max(df$abcdScore)], df$abcdScore[which.max(df$abcdScore)], col = "blue")
points(df$start[which.max(df$abcScore)], df$abcdScore[which.max(df$abcScore)], col = "red")
points(df$start[38], df$abcdScore[38], col = "purple")
# abline(v = df$start[which.max(df$dynamics)])
abline(v = center)
abline(h = 0, lty = 2)
par(mfrow = c(1,1), pch = 1)

df$start[which.max(df$abcdScore)]

## Genomics
## ABC would predict this variant
## Look at time correlation (not great for ABC selected one)
## Dynamics has some spurious correlations though
## Combining into ABCD selects the "best" candidates
## Show genomics of the connection
## This improves assignments by ABC - model for all loci

## Zoom in the size of the plot a bit

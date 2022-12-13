## Make some hic art

## Load required packages
library(plotgardener)
library(RColorBrewer)


## Square
pdf("inst/plots/hicArt.pdf", width = 6, height = 6)

pageCreate(width = 6, height = 6, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 6, height = 6, just = c("left", "top"),
         fill = "black")
plotHicSquare(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
              x = 0.5, y = 0.5, width = 6, height = 6,
              chrom = 8, chromstart = 133250000 - 500e3, chromend = 134100000 + 500e3,
              resolution = 5e03, zrange = c(0,200), norm = "SCALE",
              palette = colorRampPalette(c("black", "violet", "white")))
dev.off()


## Rectangle
pdf("inst/plots/hicArt2.pdf", width = 10, height = 4)

pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
              x = 0, y = 0, width = 10, height = 4,
              chrom = 8, chromstart = 133250000 - 500e3, chromend = 134100000 + 500e3,
              resolution = 5e03, zrange = c(0,200), norm = "SCALE",
              palette = colorRampPalette(c("black", "violet", "white")))
dev.off()

## Rectangle magma palette
pdf("inst/plots/hicArt3.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 500e3, chromend = 134100000 + 500e3,
                 resolution = 5e03, zrange = c(0,200), norm = "SCALE",
                 palette = colorRampPalette(viridis::magma(n = 9)))
dev.off()


## Rectangle magma palette
pdf("inst/plots/hicArt4.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 1500e3, chromend = 134100000 + 1500e3,
                 resolution = 5e03, zrange = c(0,100), norm = "SCALE",
                 palette = colorRampPalette(viridis::magma(n = 9)))
dev.off()


## Rectangle magma palette
pdf("inst/plots/hicArt5.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 2500e3, chromend = 134100000 + 2500e3,
                 resolution = 10e03, zrange = c(0,200), norm = "SCALE",
                 palette = colorRampPalette(viridis::magma(n = 9)))
dev.off()


## Rectangle magma palette
pdf("inst/plots/hicArt6.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 4500e3, chromend = 134100000 + 4500e3,
                 resolution = 10e03, zrange = c(0,200), norm = "SCALE",
                 palette = colorRampPalette(viridis::magma(n = 9)))
dev.off()

## Rectangle magma palette
pdf("inst/plots/hicArt7.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 6500e3, chromend = 134100000 + 6500e3,
                 resolution = 25e03, zrange = c(0,500), norm = "SCALE",
                 palette = colorRampPalette(viridis::magma(n = 9)))
dev.off()


## Rectangle night-sky palette
pdf("inst/plots/hicArt8.pdf", width = 10, height = 4)
pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
         fill = "black")
plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                 x = 0, y = 0, width = 10, height = 4,
                 chrom = 8, chromstart = 133250000 - 6500e3, chromend = 134100000 + 6500e3,
                 resolution = 25e03, zrange = c(0,500), norm = "SCALE",
                 palette = colorRampPalette(c("#120f22", "#231f3a", "#332e61",
                                              "#634690", "#c59cde", "#f5c9e1")))
dev.off()



## Rectangle a few different palettes ----------------------------------------------------

## A few different palette options
pals <- list(viridis = viridis::viridis(n = 9),
             magma = viridis::magma(n = 9),
             plasma = viridis::plasma(n = 9),
             inferno = viridis::inferno(n = 9),
             cividis = viridis::cividis(n = 9),
             mako = viridis::mako(n = 9),
             rocket = viridis::rocket(n = 9),
             turbo = viridis::turbo(n = 9),
             YlGnBu = RColorBrewer::brewer.pal(9, 'YlGnBu'))

## Start pdf
pdf("inst/plots/hicArt9.pdf", width = 10, height = 4)

## Loop through palettes
for (i in 1:length(pals)) {

  ## Page
  pageCreate(width = 10, height = 4, showGuides = FALSE, xgrid = 0, ygrid = 0)

  ## Background
  plotRect(x = 0, y = 0, width = 10, height = 4, just = c("left", "top"),
           fill = pals[[i]][1])

  ## Hi-C
  plotHicRectangle(data = "inst/extdata/hic/MOMA_THP1_WT_4320_inter.hic",
                   x = 0, y = 0, width = 10, height = 4,
                   chrom = 8,
                   chromstart = 133250000 - 6500e3,
                   chromend = 134100000 + 6500e3,
                   resolution = 25e03,
                   zrange = c(0,500),
                   norm = "SCALE",
                   palette = colorRampPalette(pals[[i]]))

  ## Labels
  plotRect(x = 0.5-0.05, y = 0.5-0.05, width = 0.75, height = 0.25, just = c('left', 'top'),
           fill = "white")
  plotText(label = names(pals)[i], x = 0.5, y = 0.5, just = c('left', 'top'),
           fontcolor = "black")

}

dev.off()

## Load required libraries
library(hictoolsr)
library(InteractionSet)

## Define promoter regions
data("rnaVstCountsHg19")
prom <- promoters(x = rnaVstCountsHg19,
                  upstream = 2000,
                  downstream = 200)

## Enhancers = RR - prom regions
data("rrCountsHg19")
enh <- subsetByOverlaps(x = rrCountsHg19,
                        ranges = prom,
                        invert = TRUE)

## Calculate interaction pairs between enh and prom within 1Mb (~20s)
epPairs <- calcInteractionPairs(gr1 = enh,
                                gr2 = prom,
                                windowSize = 1e6,
                                stackMcols = FALSE)

## TODO:
## VST normalized, combine counts (average), zscore
## scale() on rows for RR and RNA

## Cluster
## zscore plots subset by gene cluster or


## Extract counts from LIMA timepoints
# tmp <-
#   epPairs |>
#   hictoolsr::binBedpe(res = 10e3)

epPairs |>
  regions() |>
  hictoolsr::binBedpe(res = 10e3)

extractCounts()


## Save enhancer promoter table
# save(epPairs, file = "data/epPairs.rda")

## Stack metadata columns
# stackMcols(epPairs, anchor = "both")



## Define parameters
bedpe = epPairs
res = 10e3
# pos = c('center', 'start') # list('center', 1:10) | 'center'
pos = list("center", 1:length(epPairs))

epPairs$test <- "hi"

epPairs

binBedpe(bedpe = epPairs, res = 10e3, pos = pos)




hictoolsr:::binAnchor(regions(epPairs), p = 1:length(regions(epPairs)), res = 10e3)

library(data.table)
dt <- data.table(id = anchorIds(bedpe, 'second'),
                 p = pos[[2]])


r[dt$id[[1]]]
dt[, .(r[id]), by = id]
pos[[2]]

## Handle case where one or more anchors are shifted by a vector
anchors(bedpe, 'first')

## Extract regions
r <- bedpe |> regions()

## Select the unique set of regions
r[anchorIds(bedpe, 'second')]

## New binBedpe function
binBedpe <- function(bedpe, res, pos) {

  ## Case: only one position is provided
  if (length(pos) == 1) pos[[2]] <- pos[[1]]

  ## Case:

  ## Extract regions
  r <- bedpe |> regions()

  ## Bin anchors
  a1 <- hictoolsr:::binAnchor(a = r, p = pos[[1]], res = res)
  a2 <- hictoolsr:::binAnchor(a = r, p = pos[[2]], res = res)

  ## Binned GInteractions object
  gi <- GInteractions(anchor1 = anchorIds(bedpe, 'first'),
                      anchor2 = anchorIds(bedpe, 'second') + length(r),
                      regions = c(r,r))

  ## Add back metadata
  mcols(gi) <- mcols(bedpe)

  ## Return binned bedpe
  return(gi)

}


stackMcols(epPairs)

table(regions(epPairs)$node.class)


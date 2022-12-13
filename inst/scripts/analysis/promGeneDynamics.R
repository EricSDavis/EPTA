# Compare dynamics of K27 at promoters with RNA of their genes over time

## Load libraries
library(InteractionSet)
library(GenomicRanges)

## Load enh-prom table
load("data/pgPairs.rda")

## Extract activity, contact and expression matrices for each timepoint ------------------
activityRaw <-
  as.matrix(mcols(pgPairs)[, grep("anchor1.K27_m[0-9]*_r1_RAW",
                                   colnames(mcols(pgPairs)))])
activityMat <-
  as.matrix(mcols(pgPairs)[, grep("anchor1.K27_m[0-9]*_VST",
                                   colnames(mcols(pgPairs)))])
rnaRaw <-
  as.matrix(mcols(pgPairs)[, grep("anchor2.RNA_m[0-9]*_r1_RAW",
                                   colnames(mcols(pgPairs)))])
rnaMat <-
  as.matrix(mcols(pgPairs)[, grep("anchor2.RNA_m[0-9]*_VST",
                                   colnames(mcols(pgPairs)))])

## Filtering of matrices and eQTLs -------------------------------------------------------

## Filter by raw count matrices
keep <-
  rowMins(activityRaw) > 0 &
  rowMins(rnaRaw) > 0

## Apply filter
pgPairs <- pgPairs[keep,]
activityMat <- activityMat[keep,]
rnaMat <- rnaMat[keep,]

## Normalize the activity matrix ---------------------------------------------------------

## Add pseudocount, log-transform, and scale
trans <- function(x) {
  scale(log2(x + 1))
}

## Apply transformation
activityMat <- apply(activityMat, 2, trans)

## Find minimum activity and contact value
minVal <- abs(min(activityMat)) + 1

## Transformed activity and contact matrices
a <- activityMat + minVal

## Define max abc score per timepoint
aMax <- rowMaxs(a)


## Calculate a dynamics measure ----------------------------------------------------------

## Define normalization function
# sc <- function(x){(x-min(x))/(max(x)-min(x))}
sc <- \(x){rm <- rowMins(x); (x-rm)/(rowMaxs(x)-rm)} #matrix-wide

## Define average manhattan distance function
md <- \(x, y){rowMeans(abs(x-y))}

## Shift by timepoints
scores <- list()
for (i in 0:3) {
  ## Normalize for shape
  na <- sc(a)[,1:(ncol(a)-i)]
  nrna <- sc(rnaMat)[,(1+i):ncol(rnaMat)]

  ## Calculate the row-wise distance between two matrices
  scores[[i+1]] <- md(na, nrna)
}

## Find the minimum distance across shifts
d <- rowMins(do.call(cbind, scores))

## Save as a temporary promoterDistance dataset
promDyn <- d
save(promDyn, file = "data/promoterDynamics.rda")

## Show an example
i <- 4
par(mfrow = c(1,2))
plot(density(d), main = "Dynamics distribution")
abline(v = d[i], lty = 2)

plot(sc(a)[i,], type = 'b',
     main = sprintf("Min shifted distance = %s", round(d[i], 4)),
     xlab = "Timepoints",
     ylab = "Min/Max Normalized ABC/RNA values")
lines(sc(rnaMat)[i,], type = 'b', col = 'red')
legend(x = "topright",
       legend = c("ABC scores",
                  "RNA values"),
       text.col = c("black", "red"), bty = 'n')
par(mfrow = c(1,1))


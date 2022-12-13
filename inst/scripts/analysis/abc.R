## Create an impromptu abc method

## Load libraries
library(InteractionSet)
library(GenomicRanges)
library(nullranges)

## Load enh-prom table
load("data/epCounts.rda")

## Load QTLs
load("data/naive_qtls.rda")
load("data/ifng_qtls.rda")
load("data/ciff.rda")

## Extract activity, contact and expression matrices for each timepoint ------------------
activityRaw <-
  as.matrix(mcols(epCounts)[, grep("anchor1.K27_m[0-9]*_r1_RAW",
                                   colnames(mcols(epCounts)))])
activityMat <-
  as.matrix(mcols(epCounts)[, grep("anchor1.K27_m[0-9]*_VST",
                                   colnames(mcols(epCounts)))])
contactMat <-
  as.matrix(mcols(epCounts)[, grep("LIMA_.*_[0-9]*_S_0.0.0_megaMap_inter.hic",
                                   colnames(mcols(epCounts)))])
rnaRaw <-
  as.matrix(mcols(epCounts)[, grep("anchor2.RNA_m[0-9]*_r1_RAW",
                                   colnames(mcols(epCounts)))])
rnaMat <-
  as.matrix(mcols(epCounts)[, grep("anchor2.RNA_m[0-9]*_VST",
                                   colnames(mcols(epCounts)))])

## Filtering of matrices and eQTLs -------------------------------------------------------

## Filter by raw count matrices
keep <-
  rowMins(activityRaw) > 0 &
  rowMins(rnaRaw) > 0 &
  rowMins(contactMat) > 0 &
  pairdist(epCounts) <= 1e06

## Apply filter
epCounts <- epCounts[keep,]
activityMat <- activityMat[keep,]
contactMat <- contactMat[keep,]
rnaMat <- rnaMat[keep,]

## Filter for very significant eQTLS (not necessary for ciff)
naive <- naive[naive$anchor1.pvalue <= 1e-08,]
ifng <- ifng[ifng$anchor1.pvalue <= 1e-08,]

## Filter out interchromosomal & less than 1e06 distances
naive <- naive[!is.na(pairdist(naive)) & pairdist(naive) <= 1e06,]
ifng <- ifng[!is.na(pairdist(ifng)) & pairdist(ifng) <= 1e06,]
ciff <- ciff[!is.na(pairdist(ciff)) & pairdist(ciff) <= 1e06,]

## Visualize EP distance
plot(density(pairdist(naive)), col = "blue", main = "EP distance")
lines(density(pairdist(epCounts)))
lines(density(pairdist(ciff)), col = 'orange3')
legend(x = "topright",
       legend = c("All EP Pairs",
                  "eQTL-Validated EP Pairs",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "orange3"), bty = 'n')

# pdf("inst/plots/abc_epDistance.pdf", width = 6, height = 5)
plot(density(pairdist(naive)), col = "blue", main = "EP distance")
lines(density(pairdist(epCounts)))
legend(x = "topright",
       legend = c("All EP Pairs",
                  "Naive eQTL-Validated EP Pairs"),
       text.col = c("black", "blue"), bty = 'n')
# dev.off()



## Calculate ABC -------------------------------------------------------------------------

## Add pseudocount, log-transform, and scale
trans <- function(x) {
  scale(log2(x + 1))
}

## Apply transformation
activityMat <- apply(activityMat, 2, trans)
contactMat <- apply(contactMat, 2, trans)

## Find minimum activity and contact value
minVal <- abs(min(activityMat, contactMat)) + 1

## Transformed activity and contact matrices
a <- activityMat + minVal
c <- contactMat + minVal
abc <- a * c
colnames(abc) <- gsub(".*_m([0-9]*)_.*", "ABC_\\1", colnames(abc))

## Define max abc score per timepoint
abcMax <- rowMaxs(abc)

## Visualize ABC scores amongst eQTLs ----------------------------------------------------

## Visualize ABC scores
# pdf("inst/plots/abc_histogram.pdf", width = 11, height = 3)
par(mfrow = c(1, 3))
hist(a, breaks = 100)
hist(c, breaks = 100)
hist(abc, breaks=100)
par(mfrow = c(1, 1))
# dev.off()


## Subset of eQTL-validated abc scores
abcValidNaive <- abc[countOverlaps(epCounts, naive) > 0,]
abcValidIfng <- abc[countOverlaps(epCounts, ifng) > 0,]
abcValidCiff <- abc[countOverlaps(epCounts, ciff) > 0,]



## Visualize ABC scores among eQTL-validated subsets
plot(density(abc), main = "ABC Score Distributions")
lines(density(abcValidNaive), col="blue")
lines(density(abcValidIfng), col="darkgreen")
lines(density(abcValidCiff), col="orange3")
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "darkgreen", "orange3"), bty = 'n')

## Build up the plot
# pdf("inst/plots/abc_score_distributions01.pdf", width = 6, height = 5)
plot(density(abc), main = "ABC Score Distributions")
legend(x = "topright",
       legend = c("All ABC Scores"),
       text.col = c("black"), bty = 'n')
# dev.off()

# pdf("inst/plots/abc_score_distributions02.pdf", width = 6, height = 5)
  plot(density(abc), main = "ABC Score Distributions")
  lines(density(abcValidNaive), col="blue")
  legend(x = "topright",
         legend = c("All ABC Scores",
                    "Naive eQTL-Validated ABC Scores"),
         text.col = c("black", "blue"), bty = 'n')
# dev.off()

# pdf("inst/plots/abc_score_distributions03.pdf", width = 6, height = 5)
plot(density(abc), main = "ABC Score Distributions")
lines(density(abcValidNaive), col="blue")
lines(density(abcValidIfng), col="darkgreen")
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores"),
       text.col = c("black", "blue", "darkgreen"), bty = 'n')
# dev.off()


## Per timepoint
par(mfrow=c(2,4))
lapply(1:8, function(l){
  plot(density(abc[,l]), xlim = c(0, 55))
  lines(density(abcValidNaive[,l]), col = "blue", xlim = c(0, 55))
  lines(density(abcValidIfng[,l]), col = "darkgreen", xlim = c(0, 55))
  lines(density(abcValidCiff[,l]), col = "orange3", xlim = c(0, 55))
})
par(mfrow=c(1,1))


## 24 vs 0 hr difference in ABC Scores & eQTL-validated ABC Scores
plot(density(abc[,8] - abc[,1]), main = "24-0hr ABC Scores", ylim = c(0, 0.25))
lines(density(abcValidNaive[,8] - abcValidNaive[,1]), col = "blue")
lines(density(abcValidIfng[,8] - abcValidIfng[,1]), col = "darkgreen")
lines(density(abcValidCiff[,8] - abcValidCiff[,1]), col = "orange3")
abline(v = 0, lty = 2)
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "darkgreen", "orange3"), bty = 'n')

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
  nabc <- sc(abc)[,1:(ncol(abc)-i)]
  nrna <- sc(rnaMat)[,(1+i):ncol(rnaMat)]

  ## Calculate the row-wise distance between two matrices
  scores[[i+1]] <- md(nabc, nrna)
}

## Find the minimum distance across shifts
d <- 1-rowMins(do.call(cbind, scores))

## Show an example
# pdf("inst/plots/abc_dynamics_example.pdf", width = 10, height = 5)
i <- 1
par(mfrow = c(1,2))
plot(density(d), main = "Dynamics distribution")
abline(v = d[i], lty = 2)

plot(sc(abc)[i,], type = 'b',
     main = sprintf("Dynamics Value = %s", round(d[i], 4)),
     xlab = "Timepoints",
     ylab = "Min/Max Normalized ABC/RNA values")
lines(sc(rnaMat)[i,], type = 'b', col = 'red')
legend(x = "topleft",
       legend = c("ABC scores",
                  "RNA values"),
       text.col = c("black", "red"), bty = 'n')
par(mfrow = c(1,1))
# dev.off()

## Visualize dynamics distribution among all and validated eQTLs -------------------------

plot(density(d), main = "Dynamics (distance) distribution")
lines(density(d[countOverlaps(epCounts, naive) > 0]), col = 'blue')
lines(density(d[countOverlaps(epCounts, ifng) > 0]), col = 'darkgreen')
lines(density(d[countOverlaps(epCounts, ciff) > 0]), col = 'orange3')
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "blue", "darkgreen", "orange3"), bty = 'n')

## Build it up
# pdf("inst/plots/abc_dynamics_distribution01.pdf", width = 8, height = 6)
plot(density(d), main = "Dynamics distribution")
legend(x = "topleft",
       legend = c("All EP-Pairs"),
       text.col = c("black"), bty = 'n')
# dev.off()

# pdf("inst/plots/abc_dynamics_distribution02.pdf", width = 8, height = 6)
plot(density(d), main = "Dynamics distribution")
lines(density(d[countOverlaps(epCounts, naive) > 0]), col = 'blue')
legend(x = "topleft",
       legend = c("All EP-Pairs",
                  "Naive eQTL-Validated EP-Pairs"),
       text.col = c("black", "blue"), bty = 'n')
# dev.off()

# pdf("inst/plots/abc_dynamics_distribution03.pdf", width = 8, height = 6)
plot(density(d), main = "Dynamics distribution")
lines(density(d[countOverlaps(epCounts, naive) > 0]), col = 'blue')
lines(density(d[countOverlaps(epCounts, ifng) > 0]), col = 'darkgreen')
legend(x = "topleft",
       legend = c("All EP-Pairs",
                  "Naive eQTL-Validated EP-Pairs",
                  "Ifng eQTL-Validated EP-Pairs"),
       text.col = c("black", "blue", "darkgreen"), bty = 'n')
# dev.off()


## Compare to the dynamics of promoters (Katie says eQTLs are stinky - YUCK!)
load("data/promoterDynamics.rda")
plot(density(d), main = "Dynamics (distance) distribution")
lines(density(d[countOverlaps(epCounts, naive) > 0]), col = 'blue')
lines(density(d[countOverlaps(epCounts, ifng) > 0]), col = 'darkgreen')
lines(density(d[countOverlaps(epCounts, ciff) > 0]), col = 'orange3')
lines(density(promDyn), col = 'purple2')
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Ifng eQTL-Validated ABC Scores",
                  "CRISPRi-FlowFISH Validated",
                  "Promoters"),
       text.col = c("black", "blue", "darkgreen", "orange3", "purple2"), bty = 'n')

## Save ABCD results with epCounts -------------------------------------------------------

# save(epCounts, abc, d, file = "data/abcdResults.rda")


## Generate distance-matched sets of EP Pairs---------------------------------------------

## Create a separate object for holding interactions and abc/d
epABCD <- epCounts
mcols(epABCD) <- abc
epABCD$dynamics <- d

## Add EP distance to epCounts object
epABCD$distance <- pairdist(epABCD)

## Add validated annotations
epABCD$valid <- countOverlaps(epABCD, naive) > 0

## Find distance-matced EP counts
set.seed(123)
epMatched <- matchRanges(focal = epABCD[epABCD$valid],
                         pool = epABCD[!epABCD$valid],
                         covar = ~distance,
                         method = 'strat',
                         replace = FALSE)

## Assess matching quality
overview(epMatched)
plotCovariate(epMatched)
# ggplot2::ggsave("inst/plots/abc_epDistanceMatched.pdf", width = 6, height = 5)

# pdf("inst/plots/abc_distance_corrected.pdf", width = 14, height = 6)
par(mfrow = c(1, 2))
## Visualize distance-corrected ABC scores among eQTL-validated subsets
plot(density(as.matrix(mcols(epABCD[!epABCD$valid, grep("ABC", colnames(mcols(epABCD)))]))),
     main = "ABC Score Distributions")
lines(density(as.matrix(mcols(epABCD[epABCD$valid, grep("ABC", colnames(mcols(epABCD)))]))),
      col="blue")
lines(density(as.matrix(mcols(epMatched[, grep("ABC", colnames(mcols(epMatched)))]))),
      col="purple2")
legend(x = "topright",
       legend = c("All ABC Scores",
                  "Naive eQTL-Validated ABC Scores",
                  "Distance-corrected ABC Scores"),
       text.col = c("black", "blue", "purple2"), bty = 'n')

## Visualize distance-corrected dynamics scores among eQTL-validated subsets
plot(density(epABCD$dynamics[!epABCD$valid]), main = "Dynamics distribution")
lines(density(epABCD$dynamics[epABCD$valid]), col = 'blue')
lines(density(epMatched$dynamics), col = 'purple2')
legend(x = "topleft",
       legend = c("All EP-Pairs",
                  "Naive eQTL-Validated EP-Pairs",
                  "Distance-corrected EP-Pairs"),
       text.col = c("black", "blue", "purple2"), bty = 'n')
par(mfrow=c(1,1))
# dev.off()


## Visualize ABC and dynamics ------------------------------------------------------------

plot(abcMax[1:8587], d[1:8587], xlim = range(abcMax), ylim = range(d),
     main = "Bivariate ABC vs. Dynamics",
     xlab = "Max ABC Score Per Timepoint",
     ylab = "Pairwise Manhattan Distance between ABC and RNA")
points(abcMax[countOverlaps(epCounts, naive) > 0],
       d[countOverlaps(epCounts, naive) > 0],
       col = 'blue')
legend(x = "topright",
       legend = c("All EP Pairs",
                  "Naive eQTL-Validated EP Pairs"),
       text.col = c("black", "blue"), bty = 'n')

## Again for CRISPRi-FlowFISH
plot(abcMax[1:3644], d[1:3644], xlim = range(abcMax), ylim = range(d),
     main = "Bivariate ABC vs. Dynamics",
     xlab = "Max ABC Score Per Timepoint",
     ylab = "Pairwise Manhattan Distance between ABC and RNA")
points(abcMax[countOverlaps(epCounts, ciff) > 0],
       d[countOverlaps(epCounts, ciff) > 0],
       col = 'orange3')
legend(x = "topright",
       legend = c("All EP Pairs",
                  "CRISPRi-FlowFISH Validated"),
       text.col = c("black", "orange3"), bty = 'n')

abline(v = seq(min(abcMax), max(abcMax), diff(range(abcMax)) / 10))
abline(h = seq(min(d), max(d), diff(range(d)) / 10))


## Alternative bivariate visualization ---------------------------------------------------

## Compile abc and dynamics together
df <- data.frame(abc = abcMax,
                 dyn = d,
                 abcGroup = cut(abcMax, 20, include.lowest = TRUE),
                 dynGroup = cut(d, 20, include.lowest = TRUE),
                 valid = countOverlaps(epCounts, ciff) > 0)

## Create heatmap grid
grps <- expand.grid(abc = levels(df$abcGroup), dyn = levels(df$dynGroup))

## Assign groups to the data


## Subset into bins and calculate percent validated
percValid <-
  lapply(1:nrow(grps), \(i) {
    bin <- subset(df, abcGroup == grps[i,1] & dynGroup == grps[i,2])
    (sum(bin$valid) / nrow(bin)) * 100
  }) |> unlist()


## Add info back to groups and remove NaNs
grps$percValid <- percValid
grps$percValid[is.nan(grps$percValid)] <- 0


library(ggplot2)
ggplot(grps, aes(x = abc, y = dyn, fill = percValid)) +
  scale_fill_distiller(palette = "YlGnBu") +
  geom_tile()


## Rank ----------------------------------------------------------------------------------

# df <- data.frame(abc = abc[,1],
#                  valid = valid)
# df <- df[order(df$abc, decreasing = TRUE),]
#
# df$cumulative <- cumsum(df$valid)/1:nrow(df)
#
# plot(1:3000, df$cumulative[1:3000], type = 'l')
# head(df)
#
# wilcox.test(abc[,1],validABC[,1])
#
# plot(density(log2(na.omit(abcMax * d))))
# lines(density(log2(na.omit(abcMax[valid] * d[valid]))), col = 'blue')


## Old measures --------------------------------------------------------------------------

# ## Calculate row-wise correlation between abc and rna matrices
# matrixCorr <- function(A, B) {
#
#   cA <- A - rowMeans(A)
#   cB <- B - rowMeans(B)
#   sA <- sqrt(rowMeans(cA^2))
#   sB <- sqrt(rowMeans(cB^2))
#
#   return(abs(rowMeans(cA * cB) / (sA * sB)))
#
# }
# d <- matrixCorr(abc, rnaMat)
#
# ## Use distance correlation metric
# euclidean <- function(a, b) sqrt(sum((a - b)^2))
# # d <- sapply(1:nrow(abc), \(x) euclidean(scale(abc[x,]), scale(rnaMat[x,]))) ## Takes too long
#
# matrixDist <- function(A, B) {
#   sA <- scale(A)
#   sB <- scale(B)
#   dAB <- (A - B)^2
#   rAB <- rowSums(dAB)
#   return(sqrt(rAB))
# }
#
# d <- matrixDist(abc, rnaMat)
#
# d[1]
# (scale(abc[1,]) - scale(rnaMat[1,]))^2
#
# head(abc)
# head(rnaMat)
#
# dist(rbind(abc[1,], rnaMat[1,]), method = "manhattan")
# rowSums(abc[1:2,] - rnaMat[1:2,])
#
# sum(sqrt((scale(abc[1,]) - scale(rnaMat[1,]))^2))
#
# rnaMat[1,]

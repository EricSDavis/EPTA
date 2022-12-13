## Load required libraries

## Load enh-prom table
load("data/binnedEPCounts.rda", verbose = TRUE)



## Subset data by cluster
cluster <- epCounts[epCounts$anchor2.cluster == "up.mid"]

## Add epdistance
cluster$epdist <- log2(pairdist(cluster) + 1)

## Transform interaction counts
cluster$contactFreq <- log2(cluster$LIMA_THP1_WT_LPIF_omega_S_0.0.0_megaMap_inter_30.hic + 1)

## Calculate total enhancer strength
cluster$enhStr <-
  cluster |>
  mcols() |>
  {\(x) x[grep("^anchor1.*VST$", colnames(x))]}() |>
  as.matrix() |>
  rowSums()

## Subset for looped EPs
looped <- cluster[cluster$looped]
unlooped <- cluster[!cluster$looped]

## Perform matching
library(nullranges)
set.seed(123)
mgi <- matchRanges(focal = looped,
                   pool = unlooped,
                   covar = ~ enhStr, # epdist, # contactFreq +
                   method = 'strat',
                   replace = FALSE)

## Assess matching quality
overview(mgi)

library(ggplot2)
plotPropensity(mgi, type = 'lines', log = 'x', sets = c('f', 'm', 'p')) + xlim(c(0, 0.02))

library(patchwork)
plots <- lapply(covariates(mgi), plotCovariate, x=mgi, sets = c('f', 'm', 'p'))
Reduce('/', plots)


## Plot lines
library(data.table)
loopedDT = as.data.table(mcols(mgi))
loopedDT = as.data.table(mcols(looped))

## Pull data of interest
tps = c(0, .5, 1, 1.5, 2, 4, 6)
rnaDat = as.matrix(loopedDT[,grep("^anchor2.*ZSCR$", colnames(loopedDT)), with=F][,1:7])
enhDat = as.matrix(loopedDT[,grep("^anchor1.*ZSCR$", colnames(loopedDT)), with=F][,1:7])

## Make empty plot
plot(x=tps,
     y=colMeans(rnaDat),
     ylim=c(-1.5, 1.5),
     type='n',
     xlab="Hours after LPS/IFNg",
     ylab="Z-score",
     main="up.mid")
abline(h=0, col="grey")

## Add RNA standard dev poly
polygon(x=c(tps, rev(tps)),
        y=c(colMeans(rnaDat)-(colSds(rnaDat)/sqrt(nrow(rnaDat))),
            rev(colMeans(rnaDat)+(colSds(rnaDat)/sqrt(nrow(rnaDat))))),
        border=NA,
        col=alpha("grey", 0.7))

## Add ENH standard dev poly
polygon(x=c(tps, rev(tps)),
        y=c(colMeans(enhDat)-(colSds(enhDat)/sqrt(nrow(enhDat))),
            rev(colMeans(enhDat)+(colSds(enhDat)/sqrt(nrow(enhDat))))),
        border=NA,
        col=alpha("grey", 0.7))

## Add mean RNA line
lines(x=tps,
      y=colMeans(rnaDat),
      type='o', pch=19)

## Add mean ENH line
lines(x=tps,
      y=colMeans(enhDat),
      type='o', pch=19, lty=2)





## Calculate correlations with shifting
calcCorrelation <- function(hic, rna){
  # Select shifts (in time segments, aka half-hours)
  shift = seq(-10, 10, by=1)
  distance=c()
  for (i in 1:length(shift)){
    s = shift[i]
    # For positive shifts (or 0 shifts)...
    if (s >= 0){
      # ...remove shift from end of hic/add to start of rna, and calculate distance
      d=rbind(c(hic[1:(length(hic)-s)]), c(rna[(1+s):length(rna)]))
      distance[i] = dist(x = d, method="manhattan")/ncol(d)
      # For negative shifts...
    }else if (s < 0) {
      # ...add shift to start of hic/remove from end of rna, and calculate distance
      d = rbind(c(hic[(1-s):length(hic)]), c(rna[1:(length(rna)+s)]))
      distance[i] = dist(x = d, method="manhattan")/ncol(d)
    }
  }
  localMin = which(distance==min(distance[which(diff(sign(diff(distance)))==2)+1]))
  plot(x=shift*.5, y=distance, type="l", xlab="shift", ylab="manhattan distance",col="black")
  points(x=shift*.5, y=distance, pch=21, cex=1.5, col="black", bg="grey90")
  points(x=shift[localMin]*.5, y=distance[localMin], pch=21, cex=1.5, col="black", bg="black")

  return(list(localMin=shift[localMin], distance=distance))
}

# Normalize 0-1
rnaDatNorm = t(apply(rnaDat, 1, function(r) (r-min(r))/(max(r)-min(r))))
enhDatNorm = t(apply(enhDat, 1, function(r) (r-min(r))/(max(r)-min(r))))

# Extrapolate to every 30 minutes
rnaExtrap = t(apply(rnaDatNorm, 1, function(r) {
  tmp = approx(y=r, x=c(0, .5, 1, 1.5, 2, 4, 6), xout=seq(0, 6, .5))
  return(tmp$y)}
))
enhExtrap =  t(apply(enhDatNorm, 1, function(r) {
  tmp = approx(y=r, x=c(0, .5, 1, 1.5, 2, 4, 6), xout=seq(0, 6, .5))
  return(tmp$y)}
))

# Plot just the extrapolated values: should look like the original cluster
plot(x=1:13, type="n", ylim=c(0,1), xaxt='n', bty='n', las=1, xlab="", ylab="")
lines(x=c(1:13), y=colMeans(rnaExtrap, na.rm=T), col="blue", lwd=2, lty=1, type='o', pch=19)
lines(x=c(1:13), y=colMeans(enhExtrap, na.rm=T), col="red", lwd=2, lty=1, type='o', pch=19)

calcCorrelation(colMeans(rnaExtrap), colMeans(enhExtrap))

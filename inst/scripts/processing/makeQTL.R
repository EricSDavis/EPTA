## Load required packages
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(InteractionSet)


## Function to import and convert to GInteractions ---------------------------------------

assembleQTL <- function(path) {

  ## Shorten txdb name
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  ## Read in qtl
  qtl <- fread(path)

  ## Group rsid
  qtl <- qtl[, .(rsid = .(rsid)), by = eval(names(qtl)[-19])]

  ## Convert to GRanges
  gr <- GRanges(seqnames = qtl$chromosome,
                ranges = IRanges(start = qtl$position, width = 1),
                strand = "*",
                pvalue = qtl$pvalue,
                gene_id = qtl$gene_id)
  gr$rsid <- qtl$rsid

  ## Convert seqnames to UCSC format
  seqlevelsStyle(gr) <- "UCSC"

  ## Lift coordinates form hg38 to hg19
  gr <-
    system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain") |>
    import.chain() |>
    {\(x) liftOver(gr, x)}() |>
    unlist()

  ## Add seqinfo
  seqlevels(gr) <- seqlevels(txdb)
  seqinfo(gr) <- seqinfo(txdb)


  ## Match to promoters and create GInteractions ##
  ## Define promoter regions from differential rna-seq data
  prom <-
    system.file("extdata/rna/LIMA_RNA_diffGenes.tsv",
                package = "EPTA") |>
    fread() |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                             seqinfo = seqinfo(txdb)) |>
    promoters(upstream = 2000, downstream = 200)

  ## Remove excess columns from prom
  mcols(prom) <- mcols(prom)[, c("ensembl", "hgnc")]

  ## Remove version from prom ENSEMBL ids
  prom$ensembl <- gsub("\\.[0-9]*", "", prom$ensembl)

  ## Match qtl ENSEMBL ids to promoter ENSEMBL ids
  m <- match(gr$gene_id, prom$ensembl)

  ## Filter out qtls that don't have a matching promoter in LIMA
  gr <- gr[!is.na(m)]
  m  <- m[!is.na(m)]

  ## Remove ensembl from qtl (apparently no way to do this)
  mcols(gr) <- mcols(gr)[, c("pvalue", "rsid")]

  ## Construct GInteractions object linking variants to promoter qtls
  gi <- GInteractions(anchor1 = gr,
                      anchor2 = prom[m])


  return(gi)

}


## Naive ---------------------------------------------------------------------------------

naive <- assembleQTL(path = "inst/extdata/qtl/filtered.p05_Alasoo_2018_ge_macrophage_naive.all.tsv.gz")


## IFNg ----------------------------------------------------------------------------------

ifng <- assembleQTL(path = "inst/extdata/qtl/filtered.p05_Alasoo_2018_ge_macrophage_IFNg.all.tsv.gz")


## Salmonella ----------------------------------------------------------------------------

salm <- assembleQTL(path = "inst/extdata/qtl/filtered.p05_Alasoo_2018_ge_macrophage_Salmonella.all.tsv.gz")


## IFNg + Salmonella ---------------------------------------------------------------------

ifng_salm <- assembleQTL(path = "inst/extdata/qtl/filtered.p05_Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz")


## Save data -----------------------------------------------------------------------------

# save(naive, ifng, salm, ifng_salm, file = "data/qtls.rda")

## Save individually ---------------------------------------------------------------------

# load("data/qtls.rda")

# save(naive, file = "data/naive_qtls.rda")
# save(ifng, file = "data/ifng_qtls.rda")
# save(salm, file = "data/salm_qtls.rda")
# save(ifng_salm, file = "data/ifng_salm_qtls.rda")

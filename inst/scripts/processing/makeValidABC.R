library(readxl)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(InteractionSet)

## Shorten txdb name
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## Read in CRIPSR flowfish
dat <-
  readxl::read_xlsx(path = "inst/extdata/abc_datasets/41588_2019_538_MOESM3_ESM.xlsx",
                    sheet = "Supplementary Table 6a",
                    col_names = TRUE,
                    skip = 1) |>
  as.data.table()

## Convert to GRanges object
gr <- GRanges(seqnames = dat$chr,
              ranges = IRanges(start = dat$start,
                               end = dat$end),
              gene = dat$Gene,
              pvalue = dat$`Adjusted p-value`,
              signficant = dat$Significant)

## Match to promoters and create GInteractions ##
## Define promoter regions from differential rna-seq data
prom <-
  fread("inst/extdata/rna/LIMA_RNA_diffGenes.tsv") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqinfo = seqinfo(txdb)) |>
  promoters(upstream = 2000, downstream = 200)

## Match CRISPR gene ids to promoter gene ids
m <- match(gr$gene, prom$hgnc)

## Filter out qtls that don't have a matching promoter in LIMA
gr <- gr[!is.na(m)]
m  <- m[!is.na(m)]

gr

## Remove ensembl from qtl (apparently no way to do this)
mcols(gr) <- mcols(gr)[, c("pvalue", "significant")]

## Construct GInteractions object linking variants to promoter qtls
ciff <- GInteractions(anchor1 = gr,
                      anchor2 = prom[m])

## Save object
save(ciff, file = "data/ciff.rda")

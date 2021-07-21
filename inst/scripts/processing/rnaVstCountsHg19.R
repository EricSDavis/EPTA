## LIMA RNA-seq Timecourse Data for hg19

## Load required libraries
library(data.table)
library(tximeta)
library(org.Hs.eg.db)
library(DESeq2)

## Read in sample sheet as coldata
coldata <-
  fread("inst/extdata/rna/LIMA_RNA_samplesheet.txt") |>
  as.data.frame()

## Load Linked Txome
loadLinkedTxome("inst/extdata/rna/GENCODE.v19.LinkedTxome.json")

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

## Add gene symbols (use org.Hs.eg.db)
gse <- addIds(se = gse, column = "SYMBOL", gene = TRUE)

## Package as DESeqDataSet for extracting vst-normalized counts
dds <- DESeqDataSet(se = gse, design = ~Bio_Rep + Timepoint)
cts <- dds |> vst() |> assay()

## Extract ranges and add vst-normalized counts
rna <- dds |> rowRanges()
names(rna) <- NULL
mcols(rna) <- cbind(mcols(rna), cts)

## Rename and save rna table
rnaVstCountsHg19 <- rna
save(rnaVstCountsHg19, file = "data/rnaVstCountsHg19.rda")

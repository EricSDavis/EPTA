## Generate samplesheet for rnaVstCountsHg19.R
## NOTE: Expects quant folder at "inst/extdata/rna/quant"
##       These files are located on longleaf at
##       "/proj/phanstiel_lab/Data/processed/LIMA/rna/LIMA_THP1_WT_LPIF_S/proc/quant"

## Load required libraries
library(data.table)

## List directory structure as data.table
x <-
  list.dirs("inst/extdata/rna/quant", recursive = FALSE) |>
  basename() |>
  data.table()

## Split to form sample sheet
ss <- cbind(x, x[, tstrsplit(V1, "_|\\.")])

## Remove unneeded columns
ss <- ss[, c("V1", "V6", "V8")]

## Rename columns
colnames(ss) <- c("names", "Timepoint", "Bio_Rep")

## Simplify names
ss$names <- gsub(pattern = "LIMA.*_([0-9]*)_S_([0-9]).*",
                 replacement = "RNA_\\1_BR\\2",
                 x = ss$names)

## Add quant.sf file paths
ss$files <-
  list.files(path = "inst/extdata/rna/quant",
             pattern = "LIMA_RNA",
             full.names = TRUE) |>
  file.path("/quant.sf")

## Check that quant paths are correct
file.exists(ss$files) |> all()

## Write out samplesheet
fwrite(ss, file = "inst/extdata/rna/LIMA_RNA_samplesheet.txt", sep = "\t")

#' LIMA RNA-seq, VST-Normalized Count GRanges
#'
#' RNA-seq data from macrophages stimulated with LPS/IF-G for 0, 30, 60, 90, 120, 240, 360, or 1440 minutes. The dataset is a GRanges object containing the locations of genes and metadata columns with VST-normalized counts for each timepoint and biological replicate. This data object was generated with the following scripts: `inst/scripts/processing/generateLIMArnaSamplesheet.R`, and `inst/scripts/processing/rnaVstCountsHg19.R`.
#'
#' @format a GRanges object with metadata peak counts
#'
#' @docType data
#'
#' @usage data("rnaVstCountsHg19")
#'
#'
"rnaVstCountsHg19"

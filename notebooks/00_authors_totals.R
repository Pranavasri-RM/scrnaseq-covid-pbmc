# 00_authors_totals.R

# Per-cell UMI totals from the authors' published count matrices
# (Wilk et al. 2020, GSE150728; STAR + dropEst exonic/intronic matrices).
# Used in 01_qc_exploration.ipynb to compare this pipeline's kallisto|bustools counts
# with theirs and to translate their QC cutoffs onto this counting scale.

# Before running (from the repo root):
#   mkdir -p data/authors && cd data/authors
#   curl -L -o GSE150728_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150728&format=file"
#   tar -xvf GSE150728_RAW.tar
#   gunzip -kf *.rds.gz   # GEO gzipped files that R had already compressed (double compression)
#   cd ../..
#
# Run from the repo root:
#   Rscript notebooks/00_authors_totals.R
#
# Output: data/authors/<sample>_authors_totals.csv
#   columns: barcode, authors_exon, authors_intron (total UMIs per cell)

# authors_exon is what their QC cutoffs apply to: the paper's intron matrices were
# generated separately, for trajectory analysis.

library(Matrix)

authors_dir <- "data/authors"
sheet       <- read.csv("assets/samplesheet.csv")

rds_files <- list.files(authors_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) == nrow(sheet))   # one file per library in the sample sheet

for (f in rds_files) {
  # "GSM4557328_555_2_cell.counts.matrices.rds" -> "GSM4557328" -> sample name in the sample sheet
  gsm    <- sub("_.*", "", basename(f))
  sample <- sheet$sample[sheet$gsm == gsm]
  stopifnot(length(sample) == 1)              # exactly one match in the sample sheet

  # List of three genes x cells matrices: exon, intron, spanning
  m <- readRDS(f)
  stopifnot(identical(colnames(m$exon), colnames(m$intron)))   # same cells, same order

  df <- data.frame(
    barcode        = colnames(m$exon),
    authors_exon   = colSums(m$exon),
    authors_intron = colSums(m$intron)
  )
  write.csv(df, file.path(authors_dir, paste0(sample, "_authors_totals.csv")), row.names = FALSE)
  cat(sample, ":", ncol(m$exon), "cells in the authors' matrix\n")
}
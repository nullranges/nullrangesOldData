## Load H3K27Ac Peak Data

## Load required libraries
library(data.table)
library(plyranges)
library(magrittr)

## Read in merged peak counts
chipPeaks <-
  fread("inst/extdata/chip/hg19/H3K27ac/peaks/h3k27ac_hg19_counts.txt") %>%
  as_granges(keep_mcols = T) %>%
  mutate(start = start + 1)

## Simplify column names
mcols(chipPeaks) %<>%
  set_colnames(gsub('.*WT_(.*)_S_NA_([0-9]).*', '\\1_\\2', colnames(.)))

## Add seqinfo
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqinfo(chipPeaks) <- seqinfo(txdb)[seqnames(seqinfo(chipPeaks))]

## Quantile normalize peak counts
mcols(chipPeaks) %<>%
  as.matrix() %>%
  preprocessCore::normalize.quantiles() %>%
  `colnames<-`(colnames(mcols(chipPeaks)))

## Calculate fold change (PMA/NON) in chip peaks
chipPeaks$peakFC <- 
  rowMedians(as.matrix(mcols(chipPeaks)[3:4])) /
  rowMedians(as.matrix(mcols(chipPeaks)[1:2]))

## Calculate peak strength across samples
chipPeaks$peakStrength <-
  chipPeaks %>%
  plyranges::select(1:4) %>%
  mcols() %>%
  as.matrix %>%
  rowSums()

## Change data name
H3K27acPeakCountsHg19 <- chipPeaks

## Save results to a file
save(H3K27acPeakCountsHg19, file = "data/H3K27acPeakCountsHg19.rda")

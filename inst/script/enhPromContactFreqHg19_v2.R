## Load libraries
library(data.table)
library(nullrangesData)
library(plyranges)
library(InteractionSet)
library(hictoolsr)

## Load RNA-seq differential genes & H3K27Ac peak counts
data("monoMacroDEGHg19")
data("H3K27acPeakCountsHg19")

## Assemble enhancer-promoter pairs ------------------------------------------------------

## Define gene promoters
prom <- promoters(monoMacroDEGHg19, upstream = 2000, downstream = 200)

## Define enhancers (remove K27 signal that falls within a promoter)
enh <- subsetByOverlaps(H3K27acPeakCountsHg19, prom, invert = T)

## Define 2Mb regions around gene promoters
promRegions <- 
  promoters(monoMacroDEGHg19, upstream = 0, downstream = 0) %>%
  resize(., width = 2e06, fix = "start") %>%
  suppressWarnings() %>%
  trim()

## Generate enhancer-promoter pairs for enhancers that fall within 2Mb of a gene promoter
enhPromPairs <-
  findOverlaps(enh, promRegions, type = "within") %>%
  Pairs(enh, prom, hits = .)

## Convert to GInteractions object (while adding mcols)
enhPromPairs <- GInteractions(first(enhPromPairs), second(enhPromPairs))

## Calculate the distance between enh-pro pairs
enhPromPairs$epDistance <- pairdist(enhPromPairs, type="gap")


## Collect interaction frequencies between enhancer-promoter pairs -----------------------

## Apply binning function centering anchor1 and beginning at TSS for anchor2
binnedEP <- binBedpe(bedpe = enhPromPairs,
                     res = 10000,
                     a1Pos = 'center',
                     a2Pos = 2000)

## Extract loop counts for each biorep hic file (~3.5 min run time)
## (replace with local file paths)
binnedEP <- extractCounts(bedpe = binnedEP,
                          hic = c("../hic/CI_THP1_O_1_inter_30.hic",
                                  "../hic/CI_THP1_O_2_inter_30.hic",
                                  "../hic/CI_THP1_A_1_inter_30.hic",
                                  "../hic/CI_THP1_A_2_inter_30.hic"),
                          chroms = c(1:22,'X','Y'),
                          res = 10000)

## Add individual biorep contact frequencies to enhPromPairs
enhPromPairs$NON_1_hic <- mcols(binnedEP)[["CI_THP1_O_1_inter_30.hic"]]
enhPromPairs$NON_2_hic <- mcols(binnedEP)[["CI_THP1_O_2_inter_30.hic"]]
enhPromPairs$PMA_1_hic <- mcols(binnedEP)[["CI_THP1_A_1_inter_30.hic"]]
enhPromPairs$PMA_2_hic <- mcols(binnedEP)[["CI_THP1_A_2_inter_30.hic"]]

## Add total contact frequency to enhPromPairs
enhPromPairs$contactFreq <-
  mcols(enhPromPairs)[grep('hic', colnames(mcols(enhPromPairs)))] %>%
  as.matrix() %>%
  rowSums()

## Change name
enhPromContactFreqHg19_v2 <- enhPromPairs

## Save results to a file
save(enhPromContactFreqHg19_v2, file = "data/enhPromContactFreqHg19_v2.rda")

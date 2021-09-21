## Load packages
library(nullranges)
library(GenomicRanges)
library(InteractionSet)
library(plyranges)
library(magrittr)
library(data.table)

## Load helper functions
source("inst/script/util.R")

## Loading and formatting data -----------------------------------------------------------

## GM12878 loop calls
loops <-
  "inst/extdata/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz" %>%
  fread() %>%
  as_ginteractions(keep.extra.columns = FALSE) %>%
  GenomeInfoDb::`seqlevelsStyle<-`("UCSC") %>%
  binBedpe(res = 10e3, a1Pos = 'center', a2Pos = 'center') %>%
  unique()

## Adjust for 0-based start positions
regions(loops) %<>% mutate(start = start + 1)


## Load GM12878 CTCF peaks
library(AnnotationHub)
ah <- AnnotationHub()
ctcfPeaks <- ah[["AH22521"]]

## Load CTCF motifs
data("CTCF_hg19")
ctcfMotifs <- CTCF_hg19

## Copy strand information to metadata (for downstream overlaps)
ctcfMotifs$direction <- as.character(strand(ctcfMotifs))

## Subset motifs that fall within a GM12878 CTCF peak
ctcfMotifs <- 
  find_overlaps_within(ctcfMotifs, ctcfPeaks) %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY'))


## Binning the genome and generating ctcf-bound bin-pairs --------------------------------

## Tile along hg19
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
bins <-
  TxDb.Hsapiens.UCSC.hg19.knownGene %>%
  GenomeInfoDb::keepStandardChromosomes() %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY')) %>%
  seqinfo() %>%
  tileGenome(tilewidth = 10e3, cut.last.tile.in.chrom = TRUE)

## Keep bins where CTCF is bound
bins <- subsetByOverlaps(bins, ctcfMotifs)

## Calculate all ctcf-bound bin-pairs within 1Mb
binPairs <- calcPairs(gr = bins, windowSize = 1e6, mode = 'strict')


## Annotating features -------------------------------------------------------------------

## Summarize CTCF motifs by bin
bin_df <- 
  bins %>%
  group_by_overlaps(ctcfMotifs) %>%
  summarize(ctcfSignal = sum(signalValue),
            direction = direction,
            hasPlus = any(direction %in% '+'),
            hasMinus = any(direction %in% '-'),
            n_sites = n())

## Initialize empty metadata
mcols(bins) <- DataFrame(ctcfSignal = 0,
                         direction = CharacterList('*'),
                         hasPlus = FALSE,
                         hasMinus = FALSE,
                         n_sites = 0)

## Add summarized data to bins
mcols(bins[bin_df$query]) <- bin_df[-1]

## Reconstruct bin-pairs
regions(binPairs) <- bins

## Extract each anchor of binPairs
a1 <- anchors(binPairs, 'first')
a2 <- anchors(binPairs, 'second')

## Annotate looped bin-pairs
binPairs$looped[countOverlaps(binPairs, loops) == 0] <- FALSE
binPairs$looped[countOverlaps(binPairs, loops) == 1] <- TRUE

## Add ctcfSignal from both anchors
binPairs$ctcfSignal <- log2(a1$ctcfSignal + a2$ctcfSignal + 1)

## Add total number of CTCF sites from both anchors
binPairs$n_sites <- a1$n_sites + a2$n_sites
binPairs$n_sites <- cut(x = binPairs$n_sites,
                        include.lowest = TRUE,
                        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, Inf),
                        labels = c("0", "1", "2", "3", "4", "5", "6", "7+"))

## Distance
binPairs$distance <- pairdist(binPairs)

## Denote convergent CTCF sites
binPairs$convergent <- FALSE
binPairs$convergent[(a1$hasPlus & a2$hasMinus) |
                      (a2$hasPlus & a1$hasMinus)] <- TRUE

## Rename
hg19_10kb_ctcfBoundBinPairs <- binPairs

## Save results to a file
save(hg19_10kb_ctcfBoundBinPairs, file = "data/hg19_10kb_ctcfBoundBinPairs.rda", compress = "xz")

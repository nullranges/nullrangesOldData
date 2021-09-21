## Load libraries
library(nullranges)
library(magrittr)
library(data.table)
library(plyranges)
library(InteractionSet)

## Load helper functions
source("inst/script/util.R")

## Load data -----------------------------------------------------------------------------

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


## Load CTCF peaks, DNase peaks, and DNase signal tracks from AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()

ctcfPeaks <- 
  ah[["AH22521"]] %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY'), pruning.mode = 'coarse')

dnasePeaks <- 
  ah[["AH22506"]] %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY'), pruning.mode = 'coarse')

dnaseSignalTrack <- 
  ah[["AH32865"]] %>%
  rtracklayer::import.bw() %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY'), pruning.mode = 'coarse')


## Generate 10Kb bins across hg19 --------------------------------------------------------

## Tile along hg19
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
bins <-
  TxDb.Hsapiens.UCSC.hg19.knownGene %>%
  GenomeInfoDb::keepStandardChromosomes() %>%
  GenomeInfoDb::dropSeqlevels(c('chrM', 'chrY')) %>%
  seqinfo() %>%
  tileGenome(tilewidth = 10e3, cut.last.tile.in.chrom = TRUE)


## Summarize and annotate bins with CTCF, DNase, and loop data ---------------------------

## Summarize CTCF peak data by bins
ctcfPeaksByBin <- 
  bins %>%
  group_by_overlaps(ctcfPeaks) %>%
  summarize(n_ctcf_sites = n(),
            ctcfSignal = sum(signalValue))

## Summarize DNase peak data by bins
dnasePeaksByBin <- 
  bins %>%
  group_by_overlaps(dnasePeaks) %>%
  summarize(n_dnase_sites = n())

## Summarize DNase signal data by bins
dnaseSignalTrackByBin <- 
  bins %>%
  group_by_overlaps(dnaseSignalTrack) %>%
  summarize(dnaseSignal = sum(score))


## Combine bin-wise data
ctcf_dnase_df <-
  Reduce(f = function(x,y) merge(x,y, by = 'query', all = TRUE),
         x = list(ctcfPeaksByBin, dnasePeaksByBin, dnaseSignalTrackByBin),
         accumulate = FALSE) %>%
  as.data.frame()

## Replace NA values with 0
ctcf_dnase_df[is.na(ctcf_dnase_df)] <- 0
ctcf_dnase_df <- DataFrame(ctcf_dnase_df)

## Initialize empty metadata
mcols(bins) <- DataFrame(n_ctcf_sites=0, ctcfSignal=0, n_dnase_sites=0, dnaseSignal=0)

## Add summarized data to bins
mcols(bins[ctcf_dnase_df$query]) <- ctcf_dnase_df[-1]

## Log transform signals
bins$ctcfSignal <- log2(bins$ctcfSignal+1)
bins$dnaseSignal <- log2(bins$dnaseSignal+1)

## Categorize n_dnase_sites
bins$n_dnase_sites <-
  cut(bins$n_dnase_sites,
      breaks = c(0, 1, 2, 3, Inf),
      include.lowest = TRUE,
      right = FALSE, labels = c('0', '1', '2', '3+'))

## Annotated looped CTCF peaks
bins$looped <- FALSE
bins$looped[countOverlaps(bins, loops, type = 'within') >= 1] <- TRUE

## Rename
hg19_10kb_bins <- bins

## Save results to a file
save(hg19_10kb_bins, file = "data/hg19_10kb_bins.rda", compress = "xz")

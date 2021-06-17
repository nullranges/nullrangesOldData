#' nullranges H3K27ac monocyte/macrophage peak count example data
#'
#' An example dataset containing GRanges of merged H3K27ac ChIP-seq peaks in THP-1 cells either untreated (NON, monocytes) or treated with PMA (PMA, macrophages). Metadata columns contain peak counts from each condition and biological replicate (NON_1, NON_2, PMA_1, and PMA_2)  obtained from \code{bedtools multicov}. Additional metadata columns contain summary statistics calculated in "inst/script/H3K27acPeakCountsHg19.R". ChIP-seq data was aligned to hg19.
#' 
#' @format a GRanges object with metadata peak counts
#' 
#' @docType data
#' 
#' @usage data("H3K27acPeakCountsHg19")
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
#' 
"H3K27acPeakCountsHg19"

#' nullranges H3K27ac monocyte/macrophage peak count example data
#'
#' An example dataset containing GRanges of merged H3K27ac ChIP-seq peaks in THP-1 cells either untreated (NON, monocytes) or treated with PMA (PMA, macrophages). Metadata columns contain peak counts from each condition and biological replicate (NON_1, NON_2, PMA_1, and PMA_2)  obtained from \code{bedtools multicov}. Additional metadata columns contain summary statistics calculated in "inst/script/H3K27acPeakCountsHg38.R". ChIP-seq data was aligned to hg38.
#' 
#' @format a GRanges object with metadata peak counts
#' 
#' @docType data
#' 
#' @usage data("H3K27acPeakCountsHg38")
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
#' 
"H3K27acPeakCountsHg38"

#' nullranges macrophage/monocyte differential gene example data
#'
#' An example dataset containing the results of a DESeq2 analysis of THP-1 cells either untreated (NON, monocytes) or treated with PMA (PMA, macrophages). The dataset is a GRanges object containing the locations of genes and metadata columns contain the DESeq2 results. Analysis of raw data in "inst/script/monoMacroDEGHg19.R". RNA-seq was quantified with salmon-1.4.0 to hg19.
#' 
#' @format a GRanges object with metadata peak counts
#' 
#' @docType data
#' 
#' @usage data("monoMacroDEGHg19")
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
#' 
"monoMacroDEGHg19"

#' nullranges macrophage/monocyte differential gene example data
#'
#' An example dataset containing the results of a DESeq2 analysis of THP-1 cells either untreated (NON, monocytes) or treated with PMA (PMA, macrophages). The dataset is a GRanges object containing the locations of genes and metadata columns contain the DESeq2 results. Analysis of raw data in "inst/script/monoMacroDEGHg38.R". RNA-seq was quantified with salmon-1.4.0 to hg38. The Salmon quant directories are available as a gzipped tar file on Zenodo.
#' 
#' @format a GRanges object with metadata peak counts
#' 
#' @docType data
#' 
#' @usage data("monoMacroDEGHg38")
#' 
#' @source https://zenodo.org/record/4484319#.YBdM2NZOmEI
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
#' https://zenodo.org/record/4484319#.YBdM2NZOmEI
#' 
"monoMacroDEGHg38"

#' nullranges enhancer-promoter contact frequency example data
#'
#' An example dataset containing the interaction frequencies between all putative enhancer-promoter (E-P) pairs within 2Mb. Promoter regions from the monoMacroDEGHg19 dataset and enhancer regions from the H3K27acPeakCountsHg19 dataset were paired within 2Mb regions around gene promoters. E-P pairs were then assigned to the appropriate Hi-C bin at 10kb resolution, taking the center point of each enhancer and the TSS of each promoter. KR-normalized interaction frequencies between each E-P pair were extracted from bioreps of monocyte (NON_1_hic, NON_2_hic) and macrophage (PMA_1_hic, PMA_2_hic) Hi-C files. Total contact frequency and E-P distance were also calculated. Returned is a GInteractions object with all putative E-P pairs along with their contact frequencies. Analysis script is in 'inst/script/enhPromContactFreqHg19.R'.
#' 
#' 
#' @format a GInteractions object with metadata contact frequencies and E-P distance.
#' 
#' @docType data
#' 
#' @usage data("enhPromContactFreqHg19")
#' 
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
"enhPromContactFreqHg19"


#' nullranges enhancer-promoter contact frequency example data
#'
#' An example dataset containing the interaction frequencies between all putative enhancer-promoter (E-P) pairs within 2Mb. Promoter regions from the monoMacroDEGHg19 dataset and enhancer regions from the H3K27acPeakCountsHg19 dataset were paired within 2Mb regions around gene promoters. E-P pairs were then assigned to the appropriate Hi-C bin at 10kb resolution, taking the center point of each enhancer and the TSS of each promoter. KR-normalized interaction frequencies between each E-P pair were extracted from bioreps of monocyte (NON_1_hic, NON_2_hic) and macrophage (PMA_1_hic, PMA_2_hic) Hi-C files. Total contact frequency and E-P distance were also calculated. Returned is a GInteractions object with all putative E-P pairs along with their contact frequencies. Analysis script is in 'inst/script/enhPromContactFreqHg19_v2.R'.
#' 
#' The major difference between enhPromContactFreqHg19 and enhPromContactFreqHg19_v2 is the way "2Mb around the promoter" is defined. In enhPromContactFreqHg19, 2Mb is expanded around the center of each promoter - where the promoter is defined as 2kb upstream and 200bp downstream from the TSS. In enhPromContactFreqHg19_v2, 2Mb is expanded around the TSS of each gene, rather than some point in the middle of the promoter.
#' 
#' 
#' @format a GInteractions object with metadata contact frequencies and E-P distance.
#' 
#' @docType data
#' 
#' @usage data("enhPromContactFreqHg19_v2")
#' 
#' 
#' @references 
#' Phanstiel, Douglas H., Kevin Van Bortle, Damek Spacek, Gaelen T. Hess, Muhammad Saad Shamim, Ido Machol, Michael I. Love, Erez Lieberman Aiden, Michael C. Bassik, and Michael P. Snyder. 2017. “Static and Dynamic DNA Loops Form AP-1-Bound Activation Hubs during Macrophage Development.” Molecular Cell 67 (6): 1037–48.e6.
#' 
"enhPromContactFreqHg19_v2"

#' ENCODE Unified GRCh38 Exclusion List (previously "blacklist")
#'
#' The blacklist is a unification of three separate blacklists followed by extensive manual curation generated by Anshul Kundaje, Stanford https://www.encodeproject.org/files/ENCFF356LFX/
#' The three blacklists that were unified are as follows
#' 1. The autogenerated blacklist from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz . 
#' The corresponding citation is
#' Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z   
#' 2. The autogenerated peakPass blacklist from https://github.com/ewimberley/peakPass/blob/master/blacklists/hg38/peakPass60Perc_sorted.bed
#' The corresponding citation is
#' Charles E. Wimberley and Steffen Heber.Journal of Computational Biology.Feb 2020.259-268.http://doi.org/10.1089/cmb.2019.0295
#' 3. The blacklist capturing peaks in the nuclear genome that originate from reads from mitochondrial genome taken from https://github.com/caleblareau/mitoblacklist/blob/master/peaks/hg38_peaks.narrowPeak . This blacklist was generated by Caleb Lareau.
#' 
#' 
#' @format a GRanges object
#' 
#' @docType data
#' 
#' @usage data("deny")
#' 
"deny"

#' nullranges DNase hypersensitivity peaks in A549 cell example data
#'
#' An example dataset containing narrowPeak file from ENCODE. Retrieve record with 'object[["AH22505"]]' on Annotation Hub. Analysis script is in 'inst/script/DHSA549Hg38.R'.
#' 
#' @format a GRanges object with metadata score, signal value, p/q value and peak.
#' 
#' @docType data
#' 
#' @usage data("DHSA549Hg38")
#' 
"DHSA549Hg38"

#' Predicted CTCF-motif binding sites
#'
#' Genomic coordinates of predicted CTCF binding sites with motif MA0139.1 (Jaspar), in BED format, with strand orientation (directionality of binding). Human (hg19, hg38) and mouse (mm9, mm10) genomes. The binding sites were detected using the FIMO tool of the MEME suite using default settings. Extra columns include motif name (MA0139.1), score, p-value, q-value, and the motif sequence.
#' 
#' Experimental, to be submitted as an AnnotatiohHub Bioconductor package.
#' 
#' 
#' 
#' @format a GRanges object
#' 
#' @docType data
#' 
#' @usage data("CTCF_hg19")
#' 
#' 
#' @references 
#' Dozmorov MG (2021). _CTCF_. https://github.com/mdozmorov/CTCF/CTCF - R
#' package version 0.99.0, <URL: https://github.com/mdozmorov/CTCF>.
#' 
"CTCF_hg19"

#' 10Kb bins from hg19 with GM12878 metadata annotation features
#'
#' 10Kb bins were tiled across hg19 and annotated with CTCF and DNase site features from GM12878. Feature annotations for each bin include 1) the number of CTCF sites, 2) the CTCF signal strength (from peak calls), 3) the number of DNase sites, 4) the DNase signal strength (from signal tracks), and finally 5) the presence/absence of a loop to any other bin.
#'
#' 
#' @format a GRanges object with covariate metadata
#' 
#' @docType data
#' 
#' @usage data("hg19_10kb_bins")
#' 
#' 
"hg19_10kb_bins"

#' CTCF-bound 10Kb paired genomic interactions
#'
#' 10Kb bins were tiled across hg19 then subset by those which contained CTCF sites. All pairs of CTCF-bound 10Kb bins were generated and annotated with feature overlaps from GM12878. Feature annotations include 1) presence/absence of a loop between bin-pairs, 2) the total CTCF signal from both bin-pairs, 3) the number of CTCF sites from both bin-pairs, 4) the distance between bin-pairs, and finally 4) whether a convergent set of CTCF sites exists between bin-pairs.
#' 
#' 
#' @format a GInteractions object with covariate metadata
#' 
#' @docType data
#' 
#' @usage data("hg19_10kb_ctcfBoundBinPairs")
#' 
#' 
"hg19_10kb_ctcfBoundBinPairs"
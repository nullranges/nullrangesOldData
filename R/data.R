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
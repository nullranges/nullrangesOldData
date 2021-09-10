## Define helper functions for data generating scripts (from hictoolsr)

## Convert data.frames to GInteractions --------------------------------------------------

makeGInteractionsFromDataFrame <- function(df,
                                           keep.extra.columns = TRUE,
                                           starts.in.df.are.0based = FALSE) {
  
  ## Convert data.table/data.frame to DataFrame
  if ("data.frame" %in% class(df)) {
    df <- DataFrame(df)
  } else if ("DFrame" %in% class(df)) {
    df <- df
  } else {
    stop("class(df) must be either 'data.frame', 'data.table', or 'DFrame'.")
  }
  
  ## Handle improper dimensions
  if(ncol(df) < 6) {
    stop("ncol(df) must be >= 6 and start with paired interactions (i.e. chr1, start1, end1 and chr2, start2, end2).")
  }
  
  ## Split into anchors
  a1 <- df[1:3] %>% `colnames<-`(c('seqnames', 'start', 'end'))
  a2 <- df[4:6] %>% `colnames<-`(c('seqnames', 'start', 'end'))
  
  ## Convert anchors to GRanges
  a1 <- makeGRangesFromDataFrame(a1, starts.in.df.are.0based = starts.in.df.are.0based)
  a2 <- makeGRangesFromDataFrame(a2, starts.in.df.are.0based = starts.in.df.are.0based)
  
  ## Create GInteractions object
  gi <- GInteractions(a1, a2)
  
  ## Add in metadata columns
  if (keep.extra.columns & ncol(df) > 6) {
    mcols(gi) <- df[7:ncol(df)]
  }
  
  ## Return
  return(gi)
  
}

## Bin interactions by resolution --------------------------------------------------------

binBedpe <- function(bedpe, res, a1Pos, a2Pos) {
  
  if ("data.frame" %in% class(bedpe)) {
    bedpe <- try(makeGRangesFromDataFrame(bedpe))
  }
  
  ## Extract anchors
  a1 <- anchors(bedpe, type = "first")
  a2 <- anchors(bedpe, type = "second")
  
  ## Bin anchors
  a1 <- binAnchor(a = a1, p = a1Pos, res = res)
  a2 <- binAnchor(a = a2, p = a2Pos, res = res)
  
  ## Binned GInteractions object
  gi <- GInteractions(a1, a2)
  
  ## Add back metadata
  mcols(gi) <- mcols(bedpe)
  
  ## Return binned bedpe
  return(gi)
  
}

## Calculate interaction pairs -----------------------------------------------------------

calcPairs <- function(gr, windowSize, mode = 'strict') {
  
  ## Check arguments ---------------------------------------------------------------------
  
  stopifnot(isClass('GRanges', gr))
  stopifnot(isClass('integer', windowSize))
  stopifnot(length(windowSize) == 1L)
  mode <- match.arg(mode, choices = c('strict', 'reverse'))
  
  ## Begin processing --------------------------------------------------------------------
  
  ## Sort gr
  gr <- gr %>% sort()
  
  ## Define windows by windowSize
  windows <-
    gr %>%
    resize(width = windowSize, fix = 'start') %>%
    suppressWarnings() %>%
    trim()
  
  ## Group gr by windows
  ov <-
    findOverlaps(gr, windows, type = 'within') %>%
    as.data.table() %>%
    `colnames<-`(c('gr', 'windows'))
  
  ## Set up progress bar
  mx <- round(uniqueN(ov$windows)*1.01, 0)
  pb <- txtProgressBar(min = 1, max = mx, initial = NA, style = 3)
  
  ## Iterate over combinations
  ov <- ov[, {setTxtProgressBar(pb, .GRP); .(gr[1], gr[-1])}, by = windows]
  
  ## Remove NA's (last value for each chromosome)
  ov <- na.omit(ov)
  
  ## Convert coordinates to GInteractions object
  gi <- GInteractions(anchor1 = ov$V1,
                      anchor2 = ov$V2,
                      regions = gr,
                      mode = mode)
  
  ## Finish progress bar
  setTxtProgressBar(pb, value = mx)
  close(pb)
  cat('\n')
  
  return(gi)
}
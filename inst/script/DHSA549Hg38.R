library(rtracklayer)
suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
pks <- query(ah, c("DnaseSeq", "narrowPeak", "A549", "UniPk"))[[1]]
pks <- keepStandardChromosomes(pks)
seqlevels(pks) <- setdiff(seqlevels(pks), "chrM")
hg19tohg38 <- import.chain("inst/extdata/hg19ToHg38.over.chain")
peaks <- unlist(liftOver(pks, hg19tohg38))
genome(peaks) <- "hg38"
peaks <- sortSeqlevels(peaks)
peaks <- sort(peaks)
table(seqnames(peaks))

save(peaks, file = "data/DHSA549Hg38.rda")
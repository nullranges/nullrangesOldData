library(rtracklayer)
suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
pks <- query(ah, c("DnaseSeq", "narrowPeak", "A549", "UniPk"))[[1]]
pks <- keepStandardChromosomes(pks)
seqlevels(pks) <- setdiff(seqlevels(pks), "chrM")
hg19tohg38 <- import.chain("inst/extdata/hg19ToHg38.over.chain")
DHSA549Hg38 <- unlist(liftOver(pks, hg19tohg38))
genome(DHSA549Hg38) <- "hg38"
DHSA549Hg38 <- sortSeqlevels(DHSA549Hg38)
DHSA549Hg38 <- sort(DHSA549Hg38)
table(seqnames(DHSA549Hg38))

save(DHSA549Hg38, file = "data/DHSA549Hg38.rda")

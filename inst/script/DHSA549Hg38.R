library(rtracklayer)
suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
pks <- query(ah, c("DnaseSeq", "narrowPeak", "A549", "UniPk"))[[1]]
pks <- keepStandardChromosomes(pks)
seqlevels(pks) <- setdiff(seqlevels(pks), "chrM")

# hg19tohg38 <- import.chain("inst/extdata/hg19ToHg38.over.chain")
ahub.chain <- subset(ah, rdataclass == "ChainFile" & species == "Homo sapiens")
chain <- ahub.chain[ahub.chain$title == "hg19ToHg38.over.chain.gz"]
chain <- chain[[1]]
DHSA549Hg38 <- unlist(liftOver(pks, chain))
genome(DHSA549Hg38) <- "hg38"
DHSA549Hg38 <- sortSeqlevels(DHSA549Hg38)
DHSA549Hg38 <- sort(DHSA549Hg38)
table(seqnames(DHSA549Hg38))

save(DHSA549Hg38, file = "data/DHSA549Hg38.rda")

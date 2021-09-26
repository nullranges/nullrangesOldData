suppressPackageStartupMessages(library(AnnotationHub))
ah <- AnnotationHub()
pks <- query(ah, c("DnaseSeq", "narrowPeak", "A549", "UniPk"))[[1]]
pks <- keepStandardChromosomes(pks)
seqlevels(pks) <- setdiff(seqlevels(pks), "chrM")

# liftover to hg38
ah_chains <- subset(ah, rdataclass == "ChainFile" & species == "Homo sapiens")
chain <- query(ah_chains, "hg19ToHg38.over.chain.gz")[[1]]
DHSA549Hg38 <- unlist(liftOver(pks, chain))
genome(DHSA549Hg38) <- "hg38"
DHSA549Hg38 <- sortSeqlevels(DHSA549Hg38)
DHSA549Hg38 <- sort(DHSA549Hg38)
table(seqnames(DHSA549Hg38))

# add seqlengths
hg38_si <- Seqinfo(genome="hg38")
hg38_sl <- seqlengths(hg38_si)[seqlevels(DHSA549Hg38)]
seqlengths(DHSA549Hg38) <- hg38_sl

save(DHSA549Hg38, file = "data/DHSA549Hg38.rda", compress = "xz")

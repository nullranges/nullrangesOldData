library(rtracklayer)
# download from https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
exclude <- import("ENCFF356LFX.bed.gz")
genome(exclude) <- "hg38"
exclude <- sortSeqlevels(exclude)
exclude <- sort(exclude)
table(seqnames(exclude))

save(exclude, file = "data/exclude.rda")

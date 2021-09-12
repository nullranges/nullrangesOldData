library(rtracklayer)
# download from https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
exclude <- import("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/ENCFF356LFX.bed.gz")
genome(exclude) <- "hg38"
exclude <- sortSeqlevels(exclude)
exclude <- sort(exclude)
table(seqnames(exclude))

save(exclude, file = "data/exclude.rda")

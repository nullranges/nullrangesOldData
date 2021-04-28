library(rtracklayer)
# download from https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
deny <- import("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/ENCFF356LFX.bed.gz")
genome(deny) <- "hg38"
deny

save(deny, file = "data/deny.rda")

# ---- Step 1: Get data files and metadata -----
library(GEOquery)

gse <- "GSE281759"

# get metadata that contains variables + other info about samples
obj <- getGEO(gse)
obj <- obj[[1]]
meta <- pData(obj)

# save meta file so we can load it in later if we need
write.csv(meta, "GSE281759_metadata.csv")

# get main data
files <- getGEOSuppFiles(gse)

# unpack main data
utils::untar("GSE281759/GSE281759_RAW.tar") 

# ---- Step 2: Analyze data ----


BiocManager::install("GEOquery")
BiocManager::install("hgu133a.db")

library(GEOquery)
library(hgu133a.db)
library(gprofiler2)

# Parse the GEO  study of interest
gset <- getGEO("GSE2361", GSEMatrix =TRUE, getGPL=FALSE)[[1]]
dim(gset)
#Check if data are normalized
boxplot(exprs(gset),outline=FALSE)
hist(log2(exprs(gset)))
exprs(gset) <- log2(exprs(gset))
# Get the features & samples num.
dim(gset)
# Get access in all probe names.
str(gset)
# Check
summary(exprs(gset))
exprs(gset)

# Save the tissue series
tis=c(gset@phenoData@data[["description"]])


# Convert probe ids into gene symbols
geneSymbols <- mapIds(hgu133a.db, keys=rownames(exprs(gset)), column="ENSEMBL", keytype="PROBEID", multiVals="filter")
geneSymbols
# Convert the result into df
df <- stack(geneSymbols)
# Copy the gene symbols into the geo matrix
new_m <- cbind(df$values, exprs(gset))
# Make gene symbols rownames instead of the probe ids
colnames(new_m)[1] <- 'GENE'
rownames(new_m) <- new_m[,1]
new_m=new_m[,-1]
colnames(new_m) <- tis

# In order to look our list of genes in interest:
setwd("/Users/artemiskounalake")
cgo_ncgo <- read.table(file = "c_nuni.txt", sep = "\t")
cgo_ncgo <- lapply(cgo_ncgo, gsub, pattern = "\\..*", replacement = "")
v1=unlist(cgo_ncgo[1])
v1 = gconvert(v1, organism = "hsapiens", target="ENSG")
v1=unlist(v1["target"])
names(v1)=NULL
v1
cgo_ncgo1 <- read.table(file = "nc_nuni.txt", sep = "\t")
cgo_ncgo1 <- lapply(cgo_ncgo1, gsub, pattern = "\\..*", replacement = "")
v2=unlist(cgo_ncgo1[1])
v2 = gconvert(v2, organism = "hsapiens", target="ENSG")
v2=unlist(v2["target"])
names(v2)=NULL
write.table(v1, file = "c_nuni_2971.txt", row.names =FALSE, col.names = FALSE )
write.table(v2, file = "nc_nuni_2971.txt", row.names =FALSE, col.names = FALSE )
write.table(new_m, file="data_expr.txt", row.names=TRUE, col.names=TRUE)

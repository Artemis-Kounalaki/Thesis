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
# Keep for every same gene the row with max mean
Amean=rowMeans(exprs(gset))
matrix2= cbind(Amean,exprs(gset) )
matrix2=matrix2[is.element(rownames(matrix2), df$ind),]
matrix3= cbind(df$values, matrix2)
colnames(matrix3)[1] <- 'GENE'
matrix3=matrix3[order(matrix3[,1],matrix3[,2],decreasing=TRUE),]
matrix3 = matrix3[!duplicated(matrix3[,1]),]
# Make gene symbols rownames instead of the probe ids
rownames(matrix3) <- matrix3[,1]
new_m=matrix3[,-2]
new_m=new_m[,-1]
colnames(new_m) <- tis


# In order to look our list of genes in interest:
c_uni <- read.table(file = "c_uni.txt", sep = "\t")
c_uni <- lapply(c_uni, gsub, pattern = "\\..*", replacement = "")
v1=unlist(c_uni[1])
v1 = gconvert(v1, organism = "hsapiens", target="ENSG")
v1=unlist(v1["target"])
names(v1)=NULL
v1

nc_uni <- read.table(file = "nc_uni.txt", sep = "\t")
nc_uni <- lapply(nc_uni, gsub, pattern = "\\..*", replacement = "")
v2=unlist(nc_uni[1])
v2 = gconvert(v2, organism = "hsapiens", target="ENSG")
v2=unlist(v2["target"])
names(v2)=NULL

c_nuni <- read.table(file = "c_nuni.txt", sep = "\t")
c_nuni <- lapply(c_nuni, gsub, pattern = "\\..*", replacement = "")
v3=unlist(c_nuni[1])
v3 = gconvert(v3, organism = "hsapiens", target="ENSG")
v3=unlist(v3["target"])
names(v3)=NULL

nc_nuni <- read.table(file = "nc_nuni.txt", sep = "\t")
nc_nuni <- lapply(nc_nuni, gsub, pattern = "\\..*", replacement = "")
v4=unlist(nc_nuni[1])
v4 = gconvert(v4, organism = "hsapiens", target="ENSG")
v4=unlist(v4["target"])
names(v4)=NULL

write.table(v1, file = "c_uni_2971.txt", row.names =FALSE, col.names = FALSE )
write.table(v2, file = "nc_uni_2971.txt", row.names =FALSE, col.names = FALSE )
write.table(v3, file = "c_nuni_2971.txt", row.names =FALSE, col.names = FALSE )
write.table(v4, file = "nc_nuni_2971.txt", row.names =FALSE, col.names = FALSE )

write.table(new_m, file="data_expr_2971.txt", row.names=TRUE, col.names=TRUE)

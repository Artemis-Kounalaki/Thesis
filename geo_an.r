BiocManager::install("GEOquery")
library(GEOquery)

# Parse the GEO  study of interest
gset <- getGEO("GSE2361", GSEMatrix =TRUE, getGPL=FALSE)[[1]]
# Get the features & samples num.
dim(gset)
# Get access in all probe names.
str(gset)
# Check for normalization
summary(exprs(gset))
# Log transform the table
exprs(gset) <- log2(exprs(gset))
exprs(gset)

# See now the normalised counts
boxplot(exprs(gset),outline=FALSE)
# Save the tissue series
tis=c(gset@phenoData@data[["description"]])




# Install the human database with probe ids
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
# Convert probe ids into gene symbols
geneSymbols <- select(hgu133plus2.db, keys=rownames(exprs(gset)), column="SYMBOL", keytype="PROBEID", multiVals="filter")
# Convert the result into df
df <- data.frame(matrix(unlist(geneSymbols), ncol =length(geneSymbols)))
# Remove dublicates
df = df[!duplicated(df$X1),]

# Copy the gene symbols into the geo matrix
new_m <- cbind(df$X2, exprs(gset))
# Make gene symbols rownames instead of the probe ids
colnames(new_m)[1] <- 'GENE'
rownames(new_m) <- new_m[,1]
new_m=new_m[,-1]
colnames(new_m) <- tis
# In order to look our list of genes in interest:

library(gprofiler2)
# Parse the txt file as a table with proteins of interest
cons_par0 <- read.table(file = "cons_par0.txt", sep = "\n")
# Make it a vector and remove the after . numbers
cons_par0=unlist(cons_par0)
names(cons_par0)=NULL
cons_par0=gsub("\\..*","",cons_par0)
#!!
write.csv(cons_par0, file = "CONS_1.csv")
# Convert the proteins to genes
gp_cons_par0 = gconvert(cons_par0, organism = "hsapiens", target="HGNC")
# Make genes in a vector
gp_cons_par0=unlist(gp_cons_par0["target"])
names(gp_cons_par0)=NULL
gp_cons_par0[1]
new_m=new_m[!duplicated(new_m[,1]),]
new_m=new_m[new_m[,1] %in% gp_cons_par0,]
length(new_m[1,])
nrow(new_m)

# THis step
library(gprofiler2)
cgo_ncgo <- read.table(file = "cgo_ncgo_par.txt", sep = " ")
cgo_ncgo <- lapply(cgo_ncgo, gsub, pattern = "\\..*", replacement = "")
v1=unlist(cgo_ncgo[1])
v1 = gconvert(v1, organism = "hsapiens", target="HGNC")
v1=unlist(v1["target"])
names(v1)=NULL
v2=unlist(cgo_ncgo[2])
v2 = gconvert(v2, organism = "hsapiens", target="HGNC")
v2=unlist(v2["target"])
names(v2)=NULL
write.table(v1, file = "target_cgo.txt", row.names =FALSE, col.names = FALSE )
write.table(v2, file = "target_ncgo.txt", row.names =FALSE, col.names = FALSE )
write.table(new_m, file="data_expr.txt", row.names=TRUE, col.names=TRUE)

new_m=new_m[!duplicated(new_m[,1]),]
new_m1=new_m[new_m[,1] %in% v1,]
new_m2=new_m[new_m[,1] %in% v2,]
dim(new_m1)
rownames(new_m1)=c(new_m1[,1])
rownames(new_m2)=c(new_m2[,1])
new_m1=new_m1[,2:37]
new_m2=new_m2[1:267,2:37]
new_m1=`class<-`(new_m1, 'numeric')
new_m2=`class<-`(new_m2, 'numeric')
res <- t.test(new_m1, new_m2, paired = TRUE)
res
res$p.value
res$estimate
res$conf.int

mat <- cbind(v1, v2)
#write.csv(mat, file = "mat.csv")
BiocManager::install("reshape2")
library(reshape2)
library(ggplot2)
hist(new_m1[,1])
hist(new_m2[,1])

dotchart(new_m1[,1],labels=row.names(mtcars),cex=.7,
         main="Gas Milage for Car Models",
         xlab="Miles Per Gallon")
dotchart(new_m2[,1],labels=row.names(mtcars),cex=.7,
         main="Gas Milage for Car Models",
         xlab="Miles Per Gallon")

### Code illustrating how subsetting Seurat or SingleCellExperiment objects
### after dimension reduction can cause confusion
  
### When working with dimension reduction results (like embeddings), 
### it's important to be careful with subsetting.  

### Basically, if you subset the object first and 
### then recompute the embedding, you'll get a new embedding based only on the subset
### But if you simply subset the object (without recomputing), 
### you'll get the original embedding values, just for the selected cells.  

### These two approaches give different results, 
### which can easily lead to confusion if you're not clear about what you're doing.

# Simulate super simple example -------------------------------------------
library(Matrix)

n <- 2000 ##must be divisible by 2
g <- 500
deg <- 100 ##must be smaller than g
ls <- c(10^seq(-3, 2, len = deg), rep(0.01, g - deg))

## create count matrix
## first half different from second half
counts <- Matrix(t(rbind(
  sapply(ls, function(l) rpois(n/2, l)),
  sapply(rev(ls), function(l) rpois(n/2, l))
)))


# Seurat ------------------------------------------------------------------

library(Seurat)
original <- CreateSeuratObject(counts = counts)

x <- RunPCA(SCTransform(original))
x <- x[,1:(n/2)]

y <- original[,1:(n/2)]
y <- RunPCA(SCTransform(y))

### you get very different results if you subset before or after
lim <- range(sapply(list(x,y), function(z) Embeddings(z, reduction = "pca")[,1]))
plot(Embeddings(x, reduction = "pca")[,1], Embeddings(y, reduction = "pca")[,1],
     xlim = lim, ylim = lim); abline(0,1)
     
# SingleCellExperiment ----------------------------------------------------

library(SingleCellExperiment)
library(scran)
original <- SingleCellExperiment(assays = list(counts = counts))
x <- fixedPCA(logNormCounts(original), subset.row = 1:g)
x <- x[,1:(n/2)]

y <- original[,1:(n/2)]
y <- fixedPCA(logNormCounts(y), subset.row = 1:g)

### you get different results if you subset before or after

lim <- range(sapply(list(x,y), function(z) reducedDim(z, "PCA")[,1]))
plot(reducedDim(x, "PCA")[,1], reducedDim(y, "PCA")[,1], 
     xlim = lim, ylim = lim); abline(0,1)



library(MASS)
library(Matrix)
plots <- FALSE

set.seed(2024-1-12)  # For reproducibility
p <- 5000 ## number of genes

## To mimic gene networks we will create a block diagonal covariance matrix
ps <- round(p*c(0.02, 0.02, 0.01, 0.01, 0.01, rep(0.005, 15)))

mats <- lapply(ps, function(pp){
  P <- matrix(runif(pp^2, 0, 0.33), pp, pp)
  m <- t(P) %*% diag(rgamma(pp, 5, 4), pp, pp) %*% P
  return(m)
})

## inactive genes
n1 <- sum(ps)
n0 <- p - n1

# Combine blocks into a block diagonal covariance matrix
library(Matrix)
cov_matrix <- bdiag(mats)

n_samples <- 50000  # Number of samples
## generate correlated data
## each gene has same expected value across all cells:
##this implies there are NO CLUSTERS
##note: we have n1 on genes and n0 off genes, that are independent
odata <- cbind(mvrnorm(n_samples, mu = rnorm(n1, 2, 2), Sigma = cov_matrix),
               sweep(sweep(matrix(rnorm(n_samples*n0), n_samples, n0), 2, rnorm(n0, -2, 0.75), FUN = "+"), 
                     2, rgamma(n0, 1, 1), FUN = "*"))
             
## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimtail data
K <- 10 
weights <- runif(K, 0.25, 0.75); weight <- weights/sum(weights)
means <- rnorm(K, -2, 0.5)
sds <- rep(0.05, K)
z <- sample(1:K, size = n_samples/K, replace = TRUE, prob = weights)
coverage <-  rnorm(n_samples, mean = means[z], sd = sds[z])
hist(coverage,nc = 100)
## add log offset to mimic different coverate
data <- odata + coverage
## genes on the rows for Seurat
data <- t(data)
## Generate counts using Poisson variation
counts <- matrix(rpois(length(unlist(data)), exp(unlist(data))), nrow(data), ncol(data))
if(plots){ ##checks
  min(colSums(counts))
  hist(log10(rowMeans(counts)))
  hist(log10(colSums(counts)))
  hist(coverage,nc=25)
  hist(colMeans(counts==0))
}
## Leverage sparsity and add rownames and colnames
counts <- as(counts, "dgCMatrix")
gc();gc() ## Garbage collection after making smaller object
rownames(counts) <- paste0("Gene", 1:nrow(counts))
colnames(counts) <- paste0("Cell", 1:ncol(counts))

## Follow standard Suerat pipeline with defaults
library(Seurat)
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
if(plots) ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:6)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)  
seurat_obj <- RunUMAP(seurat_obj, dims = 1:6)

library(ggplot2)
p0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") 
print(p0)

## change UMAP parameters so local structure is prioritized
seurat_obj <- RunUMAP(seurat_obj, dims = 1:6, n.neighbors = 5, min.dist = 0.001)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") 
print(p1)







library(MASS)
library(Matrix)

set.seed(2024-1-12)  # For reproducibility
p <- 2000 ## number of genes

## To mimic gene networks we will create a block diagonal covariance matrix
ps <- round(p*c(0.02, 0.02, 0.01, 0.01, 0.01, rep(0.005, 15)))

mats <- lapply(ps, function(pp){
  m <- matrix(runif(pp, 0, 0.49), pp, pp)
  m <- (m + t(m))/2 ## make it symmetric
  diag(m) <- 1 ## make it positive definite
  m <- m*rgamma(1, 5, 4) ## different cell to cell variability for each network
  return(m)
})

## unactive genes
mats <- c(mats, list(diag(pmax(0.1, rgamma(p - sum(ps), 1, 1)))))

# Combine blocks into a block diagonal covariance matrix
library(Matrix)
cov_matrix <- bdiag(mats)

n_samples <- 10000  # Number of samples
# Zero mean vector: this implies there are NO CLUSTERS
mean_vector <- rnorm(ncol(cov_matrix), 0, 1)  

## generate correlated data
odata <- mvrnorm(n = n_samples, mu = mean_vector, Sigma = cov_matrix)

## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimtail data
K <- 25 
weights <- runif(K, 0.25, 075); weight <- weights/sum(weights)
means <- rnorm(K, 0, 0.5)
sds <- rep(0.075, K)
z <- sample(1:K, size = n_samples/K, replace = TRUE, prob = weights)
coverage <-  rnorm(n_samples, mean = means[z], sd = sds[z])

## add log offset to mimic different coverate
data <- odata + coverage
## genes on the rows for Seurat
data <- t(data)
## Generate counts using Poisson variation
counts <- matrix(rpois(length(unlist(data)), exp(unlist(data))), nrow(data), ncol(data))
## Leverage sparsity and add rownames and colnames
counts <- as(counts, "dgCMatrix")
rownames(counts) <- paste0("Gene", 1:nrow(counts))
colnames(counts) <- paste0("Cell", 1:ncol(counts))

## Follow standard Suerat pipeline with defaults
library(Seurat)
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)  
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
p0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("Defaults")
print(p0)

## change UMAP parameters so local structure is prioritized
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, n.neighbors = 5, min.dist = 0.001)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("n.neighbors = 5, min.dist = 0.001")
print(p1)







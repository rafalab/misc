library(MASS)
library(Matrix)
library(ggplot2)
plots <- FALSE

set.seed(2024-1-12)  # For reproducibility
p <- 5000 ## number of genes

## To mimic gene networks we will create a block diagonal covariance matrix
ps <- round(p*c(0.05, 0.04, 0.03, 0.02, 0.01, rep(0.005, 25)))

mats <- lapply(ps, function(pp){
  P <- matrix(runif(1, 0.25, .66), pp, pp)
  diag(P) <- 1
  s <- pmin(2, sqrt(1/rgamma(pp, 2, 1))) # some genes vary more than others
  m <- P*outer(s,s)
  return(m)
})

## n0 inactive genes
n1 <- sum(ps)
n0 <- p - n1

# Combine blocks into a block diagonal covariance matrix
library(Matrix)
cov_matrix <- bdiag(mats)

n_samples <- 10000  # Number of samples
## generate correlated data
## each gene has same expected value across all cells:
##this implies there are NO CLUSTERS
##note: we have n1 on genes and n0 off genes, that are independent
odata <- cbind(mvrnorm(n_samples, mu = rnorm(n1, -1, 1.5), Sigma = cov_matrix),
               sweep(
                 sweep(
                   matrix(rnorm(n_samples*n0), n_samples, n0), 2, 
                   1/rgamma(n0, 5, .2), FUN = "*"), 2, 
                 rnorm(n0, -4, 0.5), FUN = "+"))
             
## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimtail data
K <- 10
weights <- runif(K, 0.25, 0.75); weight <- weights/sum(weights)
means <- rnorm(K, 1.75, 0.75)
sds <- rep(0.075, K)
z <- sample(1:K, size = n_samples/K, replace = TRUE, prob = weights)
coverage <-  rnorm(n_samples, mean = means[z], sd = sds[z])
## add log offset to mimic different coverate
data <- odata + coverage
## genes on the rows for Seurat
data <- t(data)
## Generate counts using Poisson variation
counts <- matrix(rpois(length(unlist(data)), exp(unlist(data))), 
                 nrow(data), ncol(data))
if(plots){ ##checks
  min(colSums(counts))
  hist(log10(colSums(counts)),nc=200)
  hist(colMeans(counts==0))
  hist(coverage,nc=200)
  hist(log10(rowMeans(counts)), nc = 100)
  hist(apply(log10(counts+1), 1, sd), nc = 100)
}
## Leverage sparsity and add rownames and colnames
counts <- as(counts, "dgCMatrix")
gc();gc() ## Garbage collection after making smaller object
rownames(counts) <- paste0("Gene", 1:nrow(counts))
colnames(counts) <- paste0("Cell", 1:ncol(counts))

## Follow standard Suerat pipeline with defaults
library(Seurat)
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
if(plots) ElbowPlot(seurat_obj)
K <- 6
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:K)
seurat_obj <- FindClusters(seurat_obj)  

## change UMAP parameters so local structure is prioritized
seurat_obj <- RunUMAP(seurat_obj, dims = 1:K, n.neighbors = 15, min.dist = 0.1)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("n.neighbors = 15, min.dist = 0.1")
print(p1)

ggsave(p1, filename = "~/Desktop/umap.png", width  = 6, height = 4)

## run with defaults
seurat_obj <- RunUMAP(seurat_obj, dims = 1:K)
p0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") 
print(p0) + ggtitle("n.neighbors = 30, min.dist = 0.3")

ggsave(p0, filename = "~/Desktop/umap-0.png", width  = 6, height = 4)






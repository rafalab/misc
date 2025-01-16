library(MASS)
library(Matrix)
library(matrixStats)
library(ggplot2)

## Get a sense of data to creata realisti simulation
plots <- FALSE
  
set.seed(2024 - 1 - 12)  # For reproducibility
p <- 15000 ## number of genes

## To mimic gene networks we will create a block diagonal covariance matrix
ps <- round(p*c(0.005, 0.004, 0.003, 0.002, 0.001, rep(0.0005, 25)))
p1 <- sum(ps)
p0 <- round((p - p1)*.75)
p00 <- p - p0 - p1
mats <- lapply(ps, function(pp){
  corr <- matrix(runif(1, 0.25, .66), pp, pp) ## correlation
  diag(corr) <- 1
  s <- sqrt(1/rgamma(pp, 4, 4)) # some genes vary more than others
  return(corr*outer(s,s))
})
# Combine blocks into a block diagonal covariance matrix
cov_matrix <- bdiag(mats)

n_samples <- 10000  # Number of samples
## generate correlated data
## each gene has same expected value across all cells:
##this implies there are NO CLUSTERS
##note: we have p1 "on" genes that are correlated with other genes
## and p0 off genes, and p00 completely off genes that are independent and have no variance. Sampling variation added later.
## the genes on the rows for Seurat
odata <- rbind(
  t(mvrnorm(n_samples, mu = runif(p1, -8, -4), Sigma = cov_matrix)),
    matrix(rnorm(n_samples*p0), p0, n_samples)*sqrt(1/rgamma(p0, 3, .7)) + rnorm(p0, -12, 1.5),
    matrix(0, p00, n_samples) + rnorm(p00, -15, 1))
             
## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimental data
K <- 3
weights <- runif(K, 0.25, 0.75); weights <- weights/sum(weights)
means <- log(10^seq(3.5, 4.5, len = K))
s <- rep(0.375, K)
z <- sample(1:K, size = n_samples, replace = TRUE, prob = weights)
coverage <-  rnorm(n_samples, mean = means[z], sd = s[z])
## add log offset to mimic different coverate
data <- sweep(odata, 2, coverage, FUN = "+")
data[data > log(200)] <- log(200) ##don't let count get too big
## Generate counts using Poisson variation
counts <- matrix(rpois(length(unlist(data)), exp(unlist(data))), 
                 nrow(data), ncol(data))
## Leverage sparsity and add rownames and colnames
counts <- as(counts, "dgCMatrix")
gc();gc() ## Garbage collection after making smaller object


rownames(counts) <- paste0("Gene", 1:nrow(counts))
colnames(counts) <- paste0("Cell", 1:ncol(counts))

## Simulation ends here

## Follow standard Suerat pipeline 
library(Seurat)
options(future.globals.maxSize = 2 * 1024^3)
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
#if(plots) 
ElbowPlot(seurat_obj, ndims = 50)
D <- 10 ## "elbows are seen at around 10 and at around 30
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:D)
seurat_obj <- FindClusters(seurat_obj)  

## change UMAP parameters so local structure is prioritized
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D, n.neighbors = 15, min.dist = 0.1)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("n.neighbors = 15, min.dist = 0.1")
print(p1)

ggsave(p1, filename = "~/Desktop/umap.png", width  = 10, height = 7.5)

## run with defaults
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D)
p0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("n.neighbors = 30, min.dist = 0.3")
print(p0) 

ggsave(p0, filename = "~/Desktop/umap-0.png", width  = 10, height = 7.5)


## CHECKS
if(plots){
  ### check that simnulations mimics a real dataset
  data("pbmc_facs", package = "fastTopics")
  cov <- rowSums(pbmc_facs$counts)
  x <- pbmc_facs$counts/cov
  mu <- colMeans(x)
  sds <- colSds(as.matrix(log(x+min(x[x>0])/2)))
  
  ## Is SD for expressed genes similar to simulation
  rafalib::mypar(1,2)
  hist(sds[mu > -8], nc = 25, xlim = c(0,3))
  hist(sqrt(diag(cov_matrix)), nc = 25, xlim = c(0,3)) 
  
  ## Is distribution of mean expression for genes similar 
  rafalib::mypar(2,1)
  hist(log(mu),nc = 100,xlim = c(-18,-3))
  hist(rowMeans(odata),nc = 100,xlim = c(-18,-3))
  
  ## miniumu coverage not too small
  print(min(colSums(counts)))
  
  rafalib::mypar(1,1)
  ## % 0s by cell
  hist(colMeans(counts==0), nc =100)
  ## log mean counts match real dataset
  rafalib::mypar(2,1)
  hist(log(rowMeans(counts)), nc = 100)
  hist(log(mu), nc = 100)
  rafalib::mypar(1,1)
  qqplot(log(mu), log(rowMeans(sweep(counts, 2, colSums(counts), FUN = "/"))));abline(0,1)
  ## sd of log counts match real dataset
  rafalib::mypar(2,1)
  hist(apply(log(counts+1), 1, sd), nc = 100, xlim = c(0,4))
  hist(sds, nc = 100, xlim = c(0,4))
  rafalib::mypar(1,1)
  qqplot(sds, apply(log(counts+1), 1, sd));abline(0,1)
  ## observed coverage seem realistic?
  ## should be between 1.000 and 100,000
  hist(log10(colSums(counts)),nc=100)
}

## If you want to try another dataset
#data("pbmc_facs", package = "fastTopics")
#ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "b_cells")
#x <- t(pbmc_facs$counts[ind,])

## GLM PCA removes coverage effect
# library(fastglmpca)
# fastglmpca::set_fastglmpca_threads(11)
# o <- sample(nrow(counts), 2000)
# gpca <- fastglmpca::fit_glmpca_pois(counts[o,], 50, control = list(maxiter = 1000))
# plot(log10(colSums(counts[o,])), gpca$V[,1])




library(MASS)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(lme4)
## Get a sense of data to creata realisti simulation
plots <- FALSE
  
## Get means and SDs from real data:

data("pbmc_facs", package = "fastTopics")
ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "cd34_filtered")
x <- as.matrix(pbmc_facs$counts[ind,])
x <- x[,colSums(x)>0]
N <- rowSums(x)
obs <- 1:nrow(x)
## stats has the mean and the standard deviation of the log normal
stats <- rbind(log(colSums(x)/sum(N)), 0)
high_ind <- which(exp(stats[1,])*median(N) >= 1) ## if average of 1 count per sample we fit glmm
res <- apply(x[,high_ind], 2, function(y){
  fit <- suppressMessages(glmer(y ~ 1 + (1|obs), family = poisson(link = "log"), offset = log(N)))
  if(is.null(fit@optinfo$conv$lme4$messages)) return(c(fixef(fit), attr(VarCorr(fit)$obs, "stddev"))) else return(c(log(sum(y)/sum(N)), 0))
})
stats[,high_ind] <- res
stats[2, stats[2,]>2] <- 2 ## cap the SD

stats <- stats[,order(stats[1,], decreasing = TRUE)]
set.seed(2024 - 1 - 12)  # For reproducibility
p <- ncol(stats)

## To mimic gene networks we will create a block diagonal covariance matrix
ps <- round(p*c(0.01, 0.008, 0.006, 0.004, 0.002, rep(0.001, 25)))
starts <- head(c(1,cumsum(ps)),-1) - 1
p1 <- sum(ps)
p0 <- p - p1
mats <- lapply(seq_along(ps), function(i){
  corr <- matrix(runif(1, 0.25, .66), ps[i], ps[i]) ## correlation
  diag(corr) <- 1
  ind <- starts[i] + seq(1, ps[i])
  return(corr*outer(stats[2,ind],stats[2,ind]))
})
# Combine blocks into a block diagonal covariance matrix
cov_matrix <- bdiag(mats)

n_samples <- 10000  # Number of samples
## generate correlated data
## each gene has same expected value across all cells:
##this implies there are NO CLUSTERS
##note: we have p1 "on" genes that are correlated with other genes
## and p0 off genes, with no covariance.
odata <- rbind(
  t(mvrnorm(n_samples, mu = stats[1, 1:p1], Sigma = cov_matrix)),
    matrix(rnorm(n_samples*p0), p0, n_samples)*stats[2,(p1+1):p] + stats[1,(p1+1):p])
## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimental data
## coverage we shoot for:
K <- 3
weights <- runif(K, 0.25, 0.75); weights <- weights/sum(weights)
means <- seq(3.5, 4.5, len = K)
s <- rep(0.15, K)
z <- sample(1:K, size = n_samples, replace = TRUE, prob = weights)
lN <-  log(10^rnorm(n_samples, mean = means[z], sd = s[z]))
## add log offset to mimic different coverate
data <- sweep(odata, 2, lN, FUN = "+")

#data[data > log(200)] <- log(200) ##don't let count get too big
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
D <- 15 ## "elbows are seen at around 10 and at around 30
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
  
  ## Is SD for expressed genes similar to simulation
  rafalib::mypar(1,2)
  hist(stats[2, stats[1,] > -8], nc = 25, xlim = c(0,2))
  hist(sqrt(diag(cov_matrix)), nc = 25, xlim = c(0,2)) 
  
  ## Is distribution of mean expression for genes similar 
  rafalib::mypar(2,1)
  hist(stats[1,],nc = 100,xlim = c(-18,-3))
  hist(rowMeans(odata),nc = 100,xlim = c(-18,-3))
  
  ## miniumu coverage not too small
  print(min(colSums(counts)))
  
  rafalib::mypar(1,1)
  ## % 0s by cell
  hist(colMeans(counts==0), nc =100)
  ## log mean counts match real dataset
  rafalib::mypar(2,1)
  hist(log(rowMeans(sweep(counts, 2, colSums(counts), FUN = "/"))), nc = 100)
  hist(stats[1,], nc = 100)
  rafalib::mypar(1,1)
  qqplot(stats[1,], log(rowMeans(sweep(counts, 2, colSums(counts), FUN = "/"))));abline(0,1)
  ## sd of log counts match real dataset
  rafalib::mypar(2,1)
  N1 <- colSums(counts)
  ss1 <- apply(log(sweep(counts, 2, N1, FUN = "/") * median(N1) + 0.5), 1, sd)
  N2 <- rowSums(x)
  ss2 <- apply(log(sweep(x, 1, N2, FUN = "/") * median(N2) + 0.5), 2, sd)
  
  
  rafalib::mypar(2,1)
  hist(ss1, nc=100, xlim = c(0,1))
  hist(ss2, nc=100, xlim = c(0,1))
  
  rafalib::mypar(1,1)
  qqplot(ss1, ss2);abline(0,1)
  ## observed coverage seem realistic?
  ## should be between 1.000 and 100,000
  hist(log10(colSums(counts)),nc=100)
}

## If you want to try another dataset
#data("pbmc_facs", package = "fastTopics")
#ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "b_cells")
#x <- t(pbmc_facs$counts[ind,])

## GLM PCA removes coverage effect
library(fastglmpca)
fastglmpca::set_fastglmpca_threads(11)
o <- sample(nrow(counts), 1000)
fgpca <- fastglmpca::fit_glmpca_pois(counts[o,], 10, control = list(maxiter = 100))
plot(log10(colSums(counts[o,])), fgpca$V[,1])

library(glmpca)
gpca <- glmpca(counts[o,], 10, fam = "nb2")
plot(log10(colSums(counts[o,])), gpca$factors[,1])

sct <- Embeddings(seurat_obj, reduction = "pca")
plot(colSums(counts), sct[,2])

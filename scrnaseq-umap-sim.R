library(MASS)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(lme4)

## If TRUE diagnostic plots are made
plots <- FALSE
  
## Get means and SDs from real data.
data("pbmc_facs", package = "fastTopics")
ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "cd34_filtered")
x <- as.matrix(pbmc_facs$counts[ind,])
x <- x[,colSums(x) > 0]
n <- rowSums(x)
N <- sum(n)
l <- colSums(x)/N
phi <- colSums((x - outer(n, l))^2/outer(n, l))/(nrow(x) - 1)

## Fit Poisson-log-normal model save 
### We don't compute SD for genes with no evidence that sample var > sample mean
### To avoid computational instability we require at least one count >= 5
ind <- which(phi  > qchisq(0.95, nrow(x) - 1)/(nrow(x) - 1) & colSums(x >= 5) > 0)
### To avoid bursty genes cap SD at 2.5
sd_max <- 2.5
obs <- as.factor(1:nrow(x)) ## random effect
glmm_res <- sapply(ind, function(i){
  y <- x[,i]
  sigma0 <- sqrt(log(phi[i]))
  mu0 <- log(l[i]) - sigma0^2/2
  fit <- suppressMessages(
    glmer(y ~ 1 + (1|obs),  offset = log(n), family = poisson(link = "log"), 
          start = list(theta = sigma0, fixef = mu0),
          control = glmerControl(optimizer = "bobyqa"))
  )
  msg <- fit@optinfo$conv$lme4$messages
  if(!is.null(msg)){  ## if no convergence use sigma=0
    if(any(grepl("failed to converge", msg))){
      return(c(log(l[i]), 0))
    }
  }
  s <- sqrt(VarCorr(fit)$obs[1,1]) 
  if(s > sd_max) return(c(log(l[i]), 0))
                        
  return(c(fixef(fit), s))
})
## To help decide on sd_max value look at this plot
if(plots) hist(log2(glmm_res[2,]),nc=100) 

## for the rest of genes we assume Poisson, or sigma = 0
stats <- rbind(log(l), 0)
stats[,ind] <- glmm_res

## index for Poisson-log-normal genes
ind1 <- which(stats[2,] > 0)

## get an idea of gene-gene correlations
if(plots){
  cc <- cor(log(x[,ind1]/N*median(N) + 0.5))
  hist(cc[upper.tri(cc)],nc=100)
  boxplot(cc[upper.tri(cc)])
  heatmap(cc)
}

### START SIMULATION
set.seed(2024 - 1 - 12)  # For reproducibility

### We simulate three kinds of genes
## 1- in networks, with biological variability.
## 2- not in networks but with biological variability and
## 3- genes showing just Poisson sampling variability

## generate cov-matrix for group 1
## randomly put genes into networks
groups <- c(rep(1,25), rep(2:6, each = 10), rep(7:15, each = 5))
p1 <- length(groups)
## To mimic gene networks we will create a block diagonal covariance matrix
inds <- split(sample(ind1, p1), groups)
mats <- lapply(inds, function(i){
  corr <- matrix(runif(1, 0.25, .66), length(i), length(i)) ## correlation
  diag(corr) <- 1
  return(corr*outer(stats[2,i],stats[2,i]))
})
# Combine blocks into a block diagonal covariance matrix
cov_matrix <- bdiag(mats)

## create an index for the rest
ind0 <- setdiff(1:ncol(x), unlist(inds))
p0 <- length(ind0)

n_samples <- 10000  # Number of samples
## generate correlated data
## each gene has same expected value across all cells:
##this implies there are NO CLUSTERS
##note: we have p1 "on" genes that are correlated with other genes
## and p0 genes, with no covariance.
odata <- rbind(
  t(mvrnorm(n_samples, mu = stats[1, unlist(inds)], Sigma = cov_matrix)),
    matrix(rnorm(n_samples*p0), p0, n_samples)*stats[2,ind0] + stats[1,ind0])
## Total coverage is different from cell to cell
## Generate data from mixture model to mimic experimental data
K <- 3
weights <- runif(K, 0.25, 0.75); weights <- weights/sum(weights)
means <- seq(3.75, 4.25, len = K)
s <- rep(0.1, K)
z <- sample(1:K, size = n_samples, replace = TRUE, prob = weights)
lN <-  log(10^rnorm(n_samples, mean = means[z], sd = s[z]))
## add log offset to mimic different coverage
data <- sweep(odata, 2, lN, FUN = "+")

## Generate counts using Poisson variation
counts <- matrix(rpois(length(data), exp(unlist(data))), nrow(data), ncol(data))
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
D <- 40 ## "elbows are seen at around 10 and at around 30
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:D)
seurat_obj <- FindClusters(seurat_obj)  

## change UMAP parameters so local structure is prioritized
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D, n.neighbors = 15, min.dist = 0.1)
plot1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("n.neighbors = 15, min.dist = 0.1")
print(plot1)

ggsave(plot1, filename = "~/Desktop/umap.png", width  = 10, height = 7.5)

## run with defaults
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D)
plot0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("n.neighbors = 30, min.dist = 0.3")
print(plot0) 

ggsave(plot0, filename = "~/Desktop/umap-0.png", width  = 10, height = 7.5)


## CHECKS
if(plots){
  
  ## Is SD for expressed genes similar to simulation
  rafalib::mypar(2,1)
  hist(stats[2, ind], nc = 25, xlim = c(0,2))
  hist(sqrt(diag(cov_matrix)), nc = 25, xlim = c(0,2)) 
  
  ## Is distribution of mean expression for genes similar 
  rafalib::mypar(2,1)
  hist(stats[1,],nc = 100,xlim = c(-18,-3))
  hist(rowMeans(odata),nc = 100,xlim = c(-18,-3))
  
  ## minimum coverage not too small
  print(min(colSums(counts)))
  
  ## maximum count not too large
  print(max(counts))
  
  rafalib::mypar(1,1)
  ## % 0s by cell realistic?
  hist(colMeans(counts==0), nc =100)
  
  ## log mean counts distribution match real dataset
  rafalib::mypar(1,1)
  qqplot(stats[1,], log(rowMeans(sweep(counts, 2, colSums(counts), FUN = "/"))));abline(0,1)
  
  ## sd of log counts match real dataset
  ## keeping in mind different sample sizes
  rafalib::mypar(2,1)
  ss1 <- apply(sweep(counts, 2, colSums(counts), FUN = "/"), 1, sd)
  ss2 <- apply(sweep(x, 1, rowSums(x), FUN = "/"), 2, sd)
   
  rafalib::mypar(2,1)
  hist(log(ss1), nc=100, xlim = c(-13,-4))
  hist(log(ss2), nc=100, xlim =  c(-13,-4))
  
  rafalib::mypar(1,1)
  qqplot(log(ss1), log(ss2));abline(0,1)
  
  ## observed coverage seem realistic?
  ## should be between 1.000 and 100,000
  hist(log10(colSums(counts)),nc=100)
  
  ## do we see correlation in simulated data
  image(cor(t(odata[1:p1,])))
}

## If you want to try another dataset
#data("pbmc_facs", package = "fastTopics")
#ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "b_cells")
#x <- t(pbmc_facs$counts[ind,])
if(FALSE){
## GLM PCA removes coverage effect
o <- sample(nrow(counts), 1000)

library(fastglmpca)
fastglmpca::set_fastglmpca_threads(11)
fgpca <- fastglmpca::fit_glmpca_pois(counts[o,], 10, control = list(maxiter = 25))
plot(log10(colSums(counts[o,])), fgpca$V[,1])

library(glmpca)
gpca <- glmpca(counts[o,], 10, fam = "nb2")
plot(log10(colSums(counts[o,])), gpca$factors[,1])
cor(log10(colSums(counts[o,])), gpca$factors[,1])

sct <- Embeddings(seurat_obj, reduction = "pca")
plot(log10(colSums(counts)), sct[,1])
cor(log10(colSums(counts)), sct[,1:5])

}


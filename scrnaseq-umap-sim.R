## To simulate scRNA-Seq data with no biological groups, 
## we generate IID multinomial samples in the following way:

## Step 1: Estimate gene-specific mean and standard deviation from real data
##         - Use a real dataset to estimate the mean and standard deviation for each gene.
##         - Fit a Poisson-log-normal distribution to genes exhibiting over-dispersion.

## Step 2: Compute correlation structure
##         - Obtain the sample correlation of log-transformed counts for over-dispersed genes.

## Step 3: Generate synthetic data
##         - Generate multivariate normal data for n_samples using the estimated means, standard deviations, and correlations.

## Step 4: Handle non-over-dispersed genes
##         - Assume a fixed mean across samples for genes that do not exhibit over-dispersion.

## Step 5: Simulate sequencing coverage
##         - Generate sequencing coverage that mimics real experimental conditions.

## Step 6: Convert data to probabilities
##         - Standardize the generated data to obtain probabilities for multinomial sampling.

## Step 7: Generate final multinomial samples
##         - Generate random multinomial count data for each sample using the computed probabilities.

library(MASS)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(lme4)

## If TRUE diagnostic plots are made
plots <- FALSE

## Load real experimental data    
data("pbmc_facs", package = "fastTopics")
ind <- which(stringr::str_extract(rownames(pbmc_facs$counts), "(?<=-)[^-]+$") == "cd34_filtered")
x <- as.matrix(pbmc_facs$counts[ind,])
x <- x[,colSums(x) > 0]
n <- rowSums(x)
N <- sum(n)
l <- colSums(x)/N
phi <- colSums((x - outer(n, l))^2/outer(n, l))/(nrow(x) - 1)
## Fit Poisson-log-normal model 
### We don't compute SD for genes with no evidence that sample var > sample mean
### To avoid computational instability we require at least one count >= 5
ind <- which(phi  > qchisq(0.95, nrow(x) - 1)/(nrow(x) - 1) & colSums(x >= 5) > 0)
### To avoid super large values cap SD so that the chance of seeing values larger than max is very small
max <- log(7500) - log(10^5) ## log of largest wanted value minus log coverage
obs <- as.factor(1:nrow(x)) ## random effect
glmm_res <- sapply(ind, function(i){
  y <- x[,i]
  sigma0 <- sqrt(log(phi[i]))
  mu0 <- log(l[i]) - sigma0^2/2
  fit <- try(suppressMessages(
    glmer(y ~ 1 + (1|obs),  offset = log(n), family = poisson(link = "log"), 
          start = list(theta = sigma0, fixef = mu0),
          control = glmerControl(optimizer = "bobyqa"))
  ))
  if(class(fit) == "try-error") return(c(log(l[i]), 0))
  msg <- fit@optinfo$conv$lme4$messages
  if(!is.null(msg)){  if(any(grepl("failed to converge", msg))){
      return(c(log(l[i]), 0))
  }}
  m <- fixef(fit); s <- sqrt(VarCorr(fit)$obs[1,1])
  return(c(m, pmin(s,(max - m)/qnorm(1-1/10^5)))) ##1/10^5 chance of reaching max
})
## To help decide on sd_max value look at this plot
if (plots){
  hist(log2(glmm_res[2,]), nc = 100) 
  plot(t(glmm_res))
}

## for the rest of genes we assume Poisson, or sigma = 0
stats <- rbind(log(l), 0)
## for overdispersed genes used fitted parameters
stats[,ind] <- glmm_res

## index and size for Poisson-log-normal genes
ind1 <- which(stats[2,] > 0)
p1 <- length(ind1)
## index and size for Poisson genes
ind0 <- which(stats[2,] == 0)
p0 <- length(ind0)

## get an estimate of gene-gene correlations
cc <- cor(log(x[,ind1]/n*median(n) + 0.5))
## use estimated correlatiosn and SDs to form covariance matrix
cov_matrix <- cc*outer(stats[2,ind1], stats[2,ind1])

### START SIMULATION
set.seed(2024 - 1 - 22)  # For reproducibility
n_samples <- 25000  # Number of samples
odata <- rbind(
  t(mvrnorm(n_samples, mu = stats[1, ind1], Sigma = cov_matrix)),
  matrix(0, p0, n_samples) + stats[1,ind0])

## Generate coverage data mimicking real data
K <- 5
weights <- runif(K, 0.25, 0.75); weights <- weights/sum(weights)
means <- seq(3, 4.5, len = K)
s <- rep(0.15, K)
z <- sample(1:K, size = n_samples, replace = TRUE, prob = weights)
n <-  round(10^rnorm(n_samples, mean = means[z], sd = s[z]))
if(plots) hist(log10(n),nc=100)
ps <- sweep(exp(odata), 2, colSums(exp(odata)), FUN = "/")
counts <- sapply(1:n_samples, function(i) rmultinom(1, n[i], ps[,i]))
counts <- as(counts, "dgCMatrix")
gc();gc() ## Garbage collection after making smaller object
if(plots) print(max(counts))

rownames(counts) <- paste0("Gene", 1:nrow(counts))
colnames(counts) <- paste0("Cell", 1:ncol(counts))

## Simulation ends here
#########################

## Follow standard Seurat pipeline 
library(Seurat)
options(future.globals.maxSize = 2 * 1024^3)
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
if (plots) ElbowPlot(seurat_obj, ndims = 50)
D <- 25 
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:D)
seurat_obj <- FindClusters(seurat_obj)  

## change UMAP parameters so local structure is prioritized 
## we use the limits of what the Seurat documentation says is reasonable
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D, n.neighbors = 5, min.dist = 0.001)
plot1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("25 PCs, n.neighbors = 5, min.dist = 0.001")
print(plot1)

ggsave(plot1, filename = "~/Desktop/umap.png", width  = 10, height = 7.5)

### Find marker genes and show them on the UMAP
markers <- sapply(c(10,12,13,17), function(cl){
  de_markers <- FindMarkers(object = seurat_obj,
                            ident.1 = as.character(cl),  # Replace with the cluster number or name
                            test.use = "wilcox",   # Default test (Wilcoxon Rank Sum Test)
                            logfc.threshold = 0.25, # Log fold change threshold
                            min.pct = 0.1 )          # Minimum fraction of cells expressing the gene
  c(rownames(de_markers)[1], de_markers$p_val_adj[1])
})

plots <- lapply(markers[1,], function(gn){
  pl <- FeaturePlot(seurat_obj, features = gn,  reduction = "umap", max.cutoff = "q95", min.cutoff = "q50")
  return(pl)
})
plot2 <- gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
ggsave(plot2, filename = "~/Desktop/markers.png", width  = 10, height = 10)

## run UMAP with PC = 10 and defaults
D <- 10
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:D)
seurat_obj <- FindClusters(seurat_obj)  
seurat_obj <- RunUMAP(seurat_obj, dims = 1:D)
plot0 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") + 
  ggtitle("10 PCs, n.neighbors = 30, min.dist = 0.3")
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
  
  ## overall distribution (take a long time)
  rafalib::mypar(1,1)
  qs <- c(0, 0.5, 0.75, 0.8, 0.9, seq(0.901, 1, 0.001))
  plot(quantile(log2(x/rowSums(x)*median(rowSums(x))+0.5), qs), 
       quantile(log2(sweep(counts, 2, colSums(counts), FUN="/")*median(colSums(counts)) + 0.5), qs))
  plot(log2(0.5+quantile(x/rowSums(x), qs)), log2(0.5+quantile(sweep(counts,2,colSums(counts)), qs)))
  
}

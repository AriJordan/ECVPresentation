setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("AriHelpers.R")

# load citations data
load("AuthorCitationWeight.Rda")
W <- weightedNetwork(tot.cite) # weighted undirected graph
n <- dim(W) # number of authors
rownames(W) <- colnames(W) <- authors

dim(W)
# Plot all authors
plotNetwork(W)

# build 15-core
converg <- FALSE
old.nrow <- nrow(W)
while(!converg){
    d <- colSums(W)
    to.keep <- which(d>=15)
    if(old.nrow==length(to.keep)){
        converg <- TRUE
    }
    old.nrow <- length(to.keep)
    W <- W[to.keep,to.keep]
}
authors <- rownames(W)

dim(W)
# Plot only authors with (citations > 15)
plotNetwork(W)

# Plot only edges with (weight > 1)
g <- plotNetwork(W, remove=1)

# Plot largest singular values
partial.SVD <- irlba(W,nv=100)
plot(1:100,partial.SVD$d)


# Choose rank K with ECV
random.est <- ECV.undirected.Rank.weighted(W,40,B=30,holdout.p=0.1,soft=FALSE,fast=TRUE)
(K <- which.min(random.est$sse))
# Stability selection?, (Meinshausen & Bühlmann, 2010)


# Choose parameter tau with ECV for spectral clustering
tau.seq <- seq(0,3,by=0.1)
#system.time(tune <- EdgeCV.REG.DC.Weight(W,h.seq,K=K,B=10,holdout.p=0.1,Arash=TRUE,fast=TRUE))
#saveRDS(tune, file = "tune.Rda")
tune <- readRDS("tune.Rda")

# Pick best tau
(best.tau <- tune$gap.min.avg)

# Apply degree regularization
d <- colSums(W)
W.reg <- W + tau.seq[best.tau]*mean(d)/n
d.reg <- colSums(W.reg)

# Laplacian
L <- t(t(W.reg/sqrt(d.reg))/sqrt(d.reg))


# Spectral clustering
# Partial singular value decomposition of Laplacian
partial.SVD <- irlba(L,nv=K)

# take K leading eigenvectors, normalize
U <- partial.SVD$u
norms <- apply(U,1,function(x)sqrt(sum(x^2)))
U.norm <- U/norms

# K-means clustering
set.seed(1)
km <- kmeans(U.norm,centers=K,iter.max=500,nstart=500)


# Separate authors into the K clusters
weighted.label <- km$cluster
weighted.cluster <- list()
for(k in 1:K){
    tmp.positions <- which(weighted.label==k)
    tmp.authors <- authors[weighted.label==k]
    tmp.degrees <- d[weighted.label==k]
    tmp.index <- sort(tmp.degrees,decreasing=TRUE,index.return=TRUE)$ix
    weighted.cluster[[k]] <- cbind(tmp.authors[tmp.index],tmp.degrees[tmp.index],tmp.positions[tmp.index])
}

# Plot the K clusters
plotClusters(g, weighted.cluster, remove=1.0)

# Turn clusters into latex table
df <- rep("",K)
for(k in 1:K){
    df[k] <- paste(weighted.cluster[[k]][1:min(7, length(weighted.cluster[[k]][,1])-1),1],collapse=", ")
}
df <- data.frame(authors=df)
library(xtable)
print(xtable(df))

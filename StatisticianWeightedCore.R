setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("AriHelpers.R")

# load citations data
load("AuthorCitationWeight.Rda")
W <- (abs(tot.cite+t(tot.cite))+abs(tot.cite-t(tot.cite)))/2 # weighted graph
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

dim(W)
# Plot only authors with (citations > 15)
plotNetwork(W)

# Plot only edges with (weight > 1)
g <- plotNetwork(W, remove=1)

PSVD <- irlba(W,nv=100)

plot(1:100,PSVD$d)

# Stability selection, (Meinshausen & Peter Bühlmann, 2010)
#result <- foreach(k = 1:100, .packages='irlba')%dopar%{
#set.seed(k)
random.est <- ECV.undirected.Rank.weighted(W,40,B=30,holdout.p=0.1,soft=FALSE,fast=TRUE)
(tmp <- which.min(random.est$sse))
#}

# Count rank occurences
#SSE.K <- unlist(result)
#(occurences <- table(SSE.K))
#plot(occurences)

# Select rank
K <- tmp
# (K <- as.numeric(names(which.max(occurences))))

# Clustering with parameter
tau.seq <- seq(0,3,by=0.1)

set.seed(500)
#system.time(tune <- EdgeCV.REG.DC.Weight(W,h.seq,K=K,B=10,holdout.p=0.1,Arash=TRUE,fast=TRUE))
#saveRDS(tune, file = "tune.Rda")
tune <- readRDS("tune.Rda")

names(tune)

tune$gap.which.min
tune$gap.min.stable
tune$gap.min.avg


d <- colSums(W)
W.reg <- W + tau.seq[tune$gap.min.avg]*mean(d)/n

d.reg <- colSums(W.reg)

L0 <- t(t(W.reg/sqrt(d.reg))/sqrt(d.reg))


PSVD <- irlba(L0,nv=K)

U <- PSVD$u

dim(U)

summary(colSums(W))


norms <- apply(U,1,function(x)sqrt(sum(x^2)))

U.norm <- U/norms

set.seed(500)
km <- kmeans(U.norm,centers=K,iter.max=500,nstart=500)


weighted.label <- km$cluster
weighted.cluster <- list()
for(k in 1:K){
    print(tmp.positions <- which(weighted.label==k))
    tmp.authors <- authors[weighted.label==k]
    tmp.degrees <- d[weighted.label==k]
    tmp.index <- sort(tmp.degrees,decreasing=TRUE,index.return=TRUE)$ix
    
    weighted.cluster[[k]] <- cbind(tmp.authors[tmp.index],tmp.degrees[tmp.index],tmp.positions[tmp.index])
}

weighted.cluster

plotGroups(g, weighted.cluster, remove=1.0)

#save(km1,SSE.K,tune,file="CitationTune.Rda")
#load("CitationTune.Rda")


df <- rep("",K)
for(k in 1:K){
    df[k] <- paste(weighted.cluster[[k]][1:min(20, length(weighted.cluster[[k]][,1])-1),1],collapse=", ")
}

df <- data.frame(authors=df)

library(xtable)
print(xtable(df))

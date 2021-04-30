require(doParallel)
registerDoParallel(cores=6)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("LeiFunc-Final.R")
source("SimonFunk.R")
source("RandomHoldout.R")
source("ArashGen.R")

# load citations data
load("AuthorCitationWeight.Rda")
W <- (abs(tot.cite+t(tot.cite))+abs(tot.cite-t(tot.cite)))/2 # weighted graph
n <- dim(W) # number of authors
rownames(W) <- colnames(W) <- authors

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

PSVD <- irlba(W,nv=100)

plot(1:100,PSVD$d)

# Stability selection, (Meinhausen & Peter Bühlmann, 2010)


result <- foreach(k = 1:100, .packages='irlba')%dopar%{
    set.seed(k)
    random.est <- ECV.undirected.Rank.weighted(W,40,B=3,holdout.p=0.1,soft=FALSE,fast=TRUE)
    tmp <- which.min(random.est$sse)
}


# Count rank occurences
SSE.K <- unlist(result)
(occurences <- table(SSE.K))
plot(occurences)

# Select rank
(K <- as.numeric(names(which.max(occurences))))

(authors <- rownames(W))
length(authors)


h.seq <- seq(0,2,by=0.1)

source("RandomHoldout-FinalSP_TestWeight.R")
set.seed(500)
#system.time(tune <- EdgeCV.REG.DC.Weight(W,h.seq,K=K,B=30,holdout.p=0.1,Arash=TRUE,fast=TRUE))
#saveRDS(tune, file = "tune.rds")
tune <- readRDS("tune.rds")

names(tune)

tune$gap.which.min
tune$gap.min.stable
tune$gap.min.avg







d <- colSums(W)
W.reg <- W + h.seq[11]*mean(d)/n

d.reg <- colSums(W.reg)

L0 <- t(t(W.reg/sqrt(d.reg))/sqrt(d.reg))


PSVD <- irlba(L0,nv=K)

U <- PSVD$u

dim(U)

summary(colSums(W))


norms <- apply(U,1,function(x)sqrt(sum(x^2)))

U.norm <- U/norms

set.seed(500)
km1 <- kmeans(U.norm,centers=20,iter.max=500,nstart=500)


weighted.label1 <- km1$cluster
weighted.cluster1 <- list()
for(k in 1:K){
    tmp.authors <- authors[weighted.label1==k]
    tmp.degrees <- d[weighted.label1==k]
    tmp.index <- sort(tmp.degrees,decreasing=TRUE,index.return=TRUE)$ix
    
    weighted.cluster1[[k]] <- cbind(tmp.authors[tmp.index],tmp.degrees[tmp.index])
}

#save(km1,SSE.K,tune,file="CitationTune.Rda")

load("CitationTune.Rda")

weighted.cluster1[[1]]  ## regression and dimensionality reduction
weighted.cluster1[[2]]  ## biostatistics - medical and genomics
weighted.cluster1[[3]]  ## High dimensional statistics-covariance estimation
weighted.cluster1[[4]]  ## Bayesian statistics
weighted.cluster1[[5]] ## decision theory and wavelets
weighted.cluster1[[6]] ## semi-parametric or nonparametrics modeling
weighted.cluster1[[7]] ## Bayesian nonparametrics
weighted.cluster1[[8]] ## Biostatistics and machine learning
weighted.cluster1[[9]] ## high dimensional statistics
weighted.cluster1[[10]] ## financial statistics
weighted.cluster1[[11]] ## unclear
weighted.cluster1[[12]] ## semi-parmatric modeling
weighted.cluster1[[13]] ## high dimensional statistics
weighted.cluster1[[14]] ## high dimensional statistics
weighted.cluster1[[15]] ## multiple testing and inference
weighted.cluster1[[16]] ## Bayesian methods and machine learning
weighted.cluster1[[17]] ## high dimensional statistics
weighted.cluster1[[18]] ## Geo-statistics and bayesian methods
weighted.cluster1[[19]] ## biostatistics
weighted.cluster1[[20]] ## functional data analysis



df <- rep("",K)
for(k in 1:K){
    df[k] <- paste(weighted.cluster1[[k]][1:7,1],collapse=", ")
}

df <- data.frame(authors=df)

library(xtable)
print(xtable(df))



df <- rep("",K)
for(k in 1:K){
    df[k] <- paste(weighted.cluster1[[k]][,1],collapse=", ")
}

df <- data.frame(authors=df)

library(xtable)
print(xtable(df))






d <- colSums(W)
W.reg <- W

d.reg <- colSums(W.reg)

L0 <- t(t(W.reg/sqrt(d.reg))/sqrt(d.reg))


PSVD <- irlba(L0,nv=20,nu=20)

U <- PSVD$u

dim(U)

summary(colSums(W))


norms <- apply(U,1,function(x)sqrt(sum(x^2)))

U.norm <- U/norms

set.seed(500)
km0 <- kmeans(U.norm,centers=20,iter.max=500,nstart=500)


weighted.label0 <- km0$cluster
weighted.cluster0 <- list()
for(k in 1:20){
    tmp.authors <- authors[weighted.label0==k]
    tmp.degrees <- d[weighted.label0==k]
    tmp.index <- sort(tmp.degrees,decreasing=TRUE,index.return=TRUE)$ix
    
    weighted.cluster0[[k]] <- cbind(tmp.authors[tmp.index],tmp.degrees[tmp.index])
}



weighted.cluster0[[1]]  ## biostatistics
weighted.cluster0[[2]]  ## regression and dimensionality reduction
weighted.cluster0[[3]]  ## high dimensional statistics
weighted.cluster0[[4]]  ## high dimensional statistics
weighted.cluster0[[5]] ## Bayesian statistics
weighted.cluster0[[6]] ## decision theory and wavelets
weighted.cluster0[[7]] ## nonparametric and semiparametric modeling
weighted.cluster0[[8]] ## Biostatistics -
weighted.cluster0[[9]] ## Bayesian statistics and machine learning
weighted.cluster0[[10]] ## Bayesian nonparametrics
weighted.cluster0[[11]] ## financial statistics
weighted.cluster0[[12]] ## Geo-statistics and bayesian methods
weighted.cluster0[[13]] ## unclear
weighted.cluster0[[14]] ## high dimensional statistics
weighted.cluster0[[15]] ## Bayesian statistics and MCMC
weighted.cluster0[[16]] ## Bayesian methods and machine learning
weighted.cluster0[[17]] ## multiple testing and inference
weighted.cluster0[[18]] ## functional data analysis, non-parametric/semi-parametric modeling
weighted.cluster0[[19]] ## unclear
weighted.cluster0[[20]] ## unclear





ss <- load("StatisticianWeightedCore-Result.Rda")


## plot

dim(W)

n <- nrow(W)


A <- matrix(0,n,n)

A[as.matrix(W)>0] <- 1

diag(A) <- 0

library(igraph)

g <- graph.adjacency(A,mode="undirected")

library(RColorBrewer)

colors <- colorRampPalette(brewer.pal(9,"YlOrRd"))(8)

d <- colSums(W)

summary(d)

quantile(d,0.99)


authors <- rownames(W)

hub.index <- which(d>330)
plot.d <- d
plot.authors <- rep(NA,n)
plot.authors[hub.index] <- authors[hub.index]

V(g)$name <- plot.authors
V(g)$size <- sqrt(d)*0.5
V(g)$label.cex <- 1

plot.color <- rep(colors[1],n)
plot.color[d>quantile(d,0.5)] <- colors[2]
plot.color[d>quantile(d,0.75)] <- colors[3]
plot.color[d>quantile(d,0.9)] <- colors[4]
plot.color[d>quantile(d,0.99)] <- colors[6]



V(g)$color <- plot.color

author.layout <- layout.fruchterman.reingold(g)



pdf("AuthorWeightNet15Core.pdf",height=12,width=12)
plot(g,layout=author.layout,vertex.label=NA,edge.width=0.8)
dev.off()




pdf("AuthorWeightNet15Core_withName.pdf",height=12,width=12)
plot(g,layout=author.layout,edge.width=0.8)
dev.off()




pdf("AuthorWeightNet15Core_withName_noMargin.pdf",height=12,width=12)
plot(g,layout=author.layout,edge.width=0.8,margin=c(0,0,0,0))
dev.off()




authors[hub.index]

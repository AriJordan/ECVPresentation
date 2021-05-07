source("LeiFunc-Final.R")
source("SimonFunk.R")
source("RandomHoldout.R")
source("ArashGen.R")
source("RandomHoldout-FinalSP_TestWeight.R")

# Build weighted network
weighted.citation.network <- function(){
  load("AuthorCitationWeight.Rda")
  W <- ((abs(tot.cite+t(tot.cite))+abs(tot.cite-t(tot.cite)))/2.0)
  rownames(W) <- colnames(W) <- authors
  return(W)
}

# Build core of network
build.core <- function(W, min.citations){
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
  return(W)
}

ECV.K <- function(W, max.K, folds=10, holdout.p=0.1){
  set.seed(2)
  random.est <- ECV.undirected.Rank.weighted(W, max.K, B=folds, holdout.p=holdout.p, soft=FALSE, fast=TRUE)
  K <- which.min(random.est$sse)
  return(K)
}

ECV.tau <- function(W, try.tau, folds=10, holdout.p=0.1){
  tau.seq <- seq(0,3,by=0.1)
  #system.time(tune <- EdgeCV.REG.DC.Weight(W,h.seq,K=K,B=10,holdout.p=0.1,Arash=TRUE,fast=TRUE))
  #saveRDS(tune, file = "tune.Rda")
  tune <- readRDS("tune.Rda")
  # Pick best tau
  best.tau <- tune$gap.min.avg # avg -> Stability selection?, (Meinshausen & Bühlmann, 2010)
  return(tau.seq[best.tau])
}

regularized.spectral.clustering <- function(W, n.clusters=K, regularization=tau){
  
  # Apply degree regularization
  d <- colSums(W)
  W.reg <- W + tau*mean(d)/n
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
  return(U.norm)
}

cluster.with.kmeans <- function(M, n.clusters){
  d <- colSums(W)
  # K-means clustering
  set.seed(1)
  km <- kmeans(M,centers=K,iter.max=500,nstart=500)
  
  weighted.label <- km$cluster
  clusters <- list()
  for(k in 1:K){
    tmp.positions <- which(weighted.label==k)
    tmp.authors <- authors[weighted.label==k]
    tmp.degrees <- d[weighted.label==k]
    tmp.index <- sort(tmp.degrees,decreasing=TRUE,index.return=TRUE)$ix
    clusters[[k]] <- cbind(tmp.authors[tmp.index],tmp.degrees[tmp.index],tmp.positions[tmp.index])
  }
  return(clusters)
}

# Plot using Fruchterman-Reingold
plot.with.fr <- function(g, hide=0){
  set.seed(5) # 5
  layout.fr <- layout_with_fr(g)
  qx1 <- quantile(layout.fr[,1],0.05)
  qx2 <- quantile(layout.fr[,1],0.95)
  qy1 <- quantile(layout.fr[,2],0.05)
  qy2 <- quantile(layout.fr[,2],0.95)
  for (v in 1:length(V(g))){
    if(layout.fr[v,1] < qx1)
      layout.fr[v,1] <- layout.fr[v,1] - (layout.fr[v,1] - qx1)/2
    if(layout.fr[v,1] > qx2)
      layout.fr[v,1] <- layout.fr[v,1] - (layout.fr[v,1] - qx2)/2
    if(layout.fr[v,2] < qy1)
      layout.fr[v,2] <- layout.fr[v,2] - (layout.fr[v,2] - qy1)/2
    if(layout.fr[v,2] > qy2)
      layout.fr[v,2] <- layout.fr[v,2] - (layout.fr[v,2] - qy2)/2
  }
  
  plot(g, edge.width=sqrt(pmax(0, E(g)$weight-hide))-0.1, layout=layout.fr, margin=c(-0.2,-1,-0.25,-1))
}

# Plot a weighted network in red-yellow
plot.network <- function(W, hide=0){
  n <- nrow(W)
  
  #A[as.matrix(W)>0] <- 1
  A <- as.matrix(W)
  diag(A) <- 0
  
  library(igraph)
  g <- graph.adjacency(A,mode="undirected")
  
  E(g)$weight <- 1
  g <- simplify(g)
  
  d <- colSums(W)
  plot.d <- d
  
  authors <- rownames(W)
  hub.index <- which(d>225)
  plot.authors <- rep(NA,n)
  plot.authors[hub.index] <- authors[hub.index]
  
  V(g)$name <- plot.authors
  V(g)$size <- sqrt(d)*0.4
  V(g)$label.cex <- 1
  
  library(RColorBrewer)
  colors <- colorRampPalette(brewer.pal(9,"YlOrRd"))(8)
  plot.color <- rep(colors[1],n)
  plot.color[d>quantile(d,0.5)] <- colors[2]
  plot.color[d>quantile(d,0.75)] <- colors[3]
  plot.color[d>quantile(d,0.9)] <- colors[4]
  plot.color[d>quantile(d,0.95)] <- colors[5]
  plot.color[d>quantile(d,0.98)] <- colors[6]
  V(g)$color <- plot.color
  
  plot.with.fr(g, hide)
  return(g)
}

plot.colored.clusters <- function(g, clusters, hide=0){
  K <- length(clusters)
  library(Polychrome)
  set.seed(2) # 3
  colors = createPalette(K,  c("#ff0000", "#ffff00", "#0000ff"))
  for(k in 1:K){
    for (a in 1:length(clusters[[k]][,3])){
      V(g)$color[as.numeric(clusters[[k]][a, 3])] <- colors[k]
    }
  }
  plot.with.fr(g, hide)
}

print.as.table <- function(clusters, authors.per.cluster){
  K <- length(clusters)
  df <- rep("",K)
  for(k in 1:K){
    df[k] <- paste(clusters[[k]][1:min(authors.per.cluster, length(clusters[[k]][,1])-1),1],collapse=", ")
  }
  df <- data.frame(authors=df)
  library(xtable)
  print(xtable(df))
}


#pdf("AuthorWeightNet15Core.pdf",height=12,width=12)

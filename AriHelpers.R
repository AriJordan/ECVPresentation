source("LeiFunc-Final.R")
source("SimonFunk.R")
source("RandomHoldout.R")
source("ArashGen.R")
source("RandomHoldout-FinalSP_TestWeight.R")

# Build weighted network
weightedNetwork <- function(tot.cite){
  return ((abs(tot.cite+t(tot.cite))+abs(tot.cite-t(tot.cite)))/2.0)
}

# Plot using Fruchterman-Reingold
plot.with.fr <- function(g, remove=0){
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
  
  plot(g, edge.width=sqrt(pmax(0, E(g)$weight-remove))-0.1, layout=layout.fr, margin=c(-0.2,-1,-0.25,-1))
}

# Plot a weighted network in red-yellow
plotNetwork <- function(W, remove=0){
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
  
  plot.with.fr(g, remove)
  return(g)
}

plotClusters <- function(g, weighted.cluster, remove=0){
  K <- length(weighted.cluster)
  library(Polychrome)
  set.seed(2) # 3
  colors = createPalette(K,  c("#ff0000", "#ffff00", "#0000ff"))
  for(k in 1:K){
    for (a in 1:length(weighted.cluster[[k]][,3])){
      V(g)$color[as.numeric(weighted.cluster[[k]][a, 3])] <- colors[k]
    }
  }
  plot.with.fr(g, remove)
}




#pdf("AuthorWeightNet15Core.pdf",height=12,width=12)

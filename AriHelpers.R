
# Plot a weighted network in red-yellow
plotNetwork <- function(W, remove=0){
  set.seed(5) # 5
  dim(W)
  
  n <- nrow(W)
  
  
  #A <- matrix(0,n,n)
  
  #A[as.matrix(W)>0] <- 1
  A <- as.matrix(W)
  diag(A) <- 0
  
  library(igraph)
  
  g <- graph.adjacency(A,mode="undirected")
  E(g)$weight <- 1
  g <- simplify(g)
  
  library(RColorBrewer)
  
  colors <- colorRampPalette(brewer.pal(9,"YlOrRd"))(8)
  
  d <- colSums(W)
  
  summary(d)
  
  quantile(d,0.99)
  
  
  authors <- rownames(W)
  
  hub.index <- which(d>225)
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
  plot.color[d>quantile(d,0.95)] <- colors[5]
  plot.color[d>quantile(d,0.98)] <- colors[6]
  
  
  
  V(g)$color <- plot.color
  
  plot(g, edge.width=sqrt(pmax(0, E(g)$weight-remove))-0.1, layout=layout_with_fr, margin=c(-0.2,-1,-0.2,-1))
  
  #author.layout <- layout.fruchterman.reingold(g)
  #plot(g,layout=author.layout,vertex.label=NA,edge.width=0.8)
}

#pdf("AuthorWeightNet15Core.pdf",height=12,width=12)

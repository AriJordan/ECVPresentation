setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("AriHelpers.R")

# Load citations data
W <- weighted.citation.network() # weighted undirected graph

# Total number of authors
length(authors <- rownames(W))

# Plot all authors
plot.network(W)

# build 15-core
W <- build.core(W, min.citations=15)

# Number of authors with (citations >= 15)
(length(authors <- rownames(W)))

# Plot only authors with (citations >= 15)
plot.network(W)

# Plot only edges with (weight > 1)
g <- plot.network(W, hide=1)

# Plot 100 largest singular values
plot(irlba(W, nv=100)$d, main= "100 biggest singular values of W", xlab="Index", ylab="Value")

# Choose rank K with ECV
(K <- ECV.K(W, max.K=40, folds=10, holdout.p=0.1))

# Choose parameter tau with ECV
(tau <- ECV.tau(W, K, try.tau=seq(0, 3, by=0.1), folds=10, holdout.p=0.1))

# Obtain clusters based on W, K and tau
clusters <- regularized.spectral.clustering.with.kmeans(W, n.clusters=K, regularization=tau)

# Plot the K clusters in colors
plot.colored.clusters(g, clusters, hide=1.0)

# Turn clusters into latex table
print.as.table(clusters, authors.per.cluster=7)

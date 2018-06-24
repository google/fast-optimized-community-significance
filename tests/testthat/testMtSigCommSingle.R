context("Overall SOCS Check")
library(igraph)
RUN_SEED <- 123454321

test_that("A random community is scored appropriately", {
  set.seed(RUN_SEED)
  network <- degree.sequence.game(rep(5, 20), method="simple")
  membership <- cluster_louvain(network)$membership
  comm <- which(membership == 1)
  adj <- get.adjacency(network)
  d <- degree(network)
  N <- length(V(network))
  edges <- get.edgelist(network)
  m <- nrow(edges)
  pval <- .SOCS(comm, adj, edges, d, N, m)
  expect_equal(signif(pval, 7), 0.3274649)
})

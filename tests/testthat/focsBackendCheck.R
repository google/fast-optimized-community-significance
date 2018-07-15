context("FOCS Backend Check")
library(igraph)
RUN_SEED <- 123454321

test_that("A random community is scored appropriately", {
  set.seed(RUN_SEED)
  network <- degree.sequence.game(rep(5, 20), method="simple")
  membership <- cluster_louvain(network)$membership
  comm <- which(membership == 1)
  score <- unlist(FOCS(list(comm), network, p=0.25))
  expect_equal(signif(score, 7), 0.2162592)
})

# Copyright 2018 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

context("FOCS Backend Check - Bipartite")
library(igraph)
RUN_SEED <- 123454321

AddClique <- function (k1, k2, edges=NULL) {

  if (is.null(edges)) {
    N <- 0
  } else {
    N <- max(edges)
  }

  new.edges <- matrix(0, ncol=2, nrow=k1 * k2)
  counter <- 1
  for (i in 1:k1) {
    for (j in (k1 + 1):(k1 + k2)) {
      new.edges[counter, ] <- c(i + N, j + N)
      counter <- counter + 1
    }
  }

  if (is.null(edges)) {
    return(new.edges)
  } else {
    return(rbind(edges, new.edges))
  }
}

ComputeGraphStatistics <- function (edges) {
  N <- max(edges[, c(1:2)])
  if (ncol(edges) == 2) {
    values <- rep(1, nrow(edges))
  } else {
    values <- edges[, 3]
  }
  adj <- Matrix::sparseMatrix(i=edges[, 1], j=edges[, 2],
                              x=values, dims=c(N, N), symmetric=TRUE)
  d <- Matrix::colSums(adj)
  m <- sum(d) / 2
  return(list(N=N, adj=adj, d=d, m=m))
}

test_that("p-score function works on small clique graph.", {
  # Compute p-scores on a 2x2 bipartite clique next to a 3x3 bipartite clique.
  # Some constant values:
  #   C' has 3 * 3 = 9 stubs.
  #   C  has 0 outgoing stubs.
  #   Removing any node adds 2 outgoing stubs to C.
  #   Every node also has 2 outgoing stubs.
  #   So the p-scores should be identically:
  #     phyper(2 - 1, 2, 9, 2, lower.tail=FALSE) = 0.01818182
  #   ...with lower bounds of zero.
  # The p-scores are tied, so worst = 7.
  toy.graph <- AddClique(2, 2, AddClique(3, 3))
  graph.stats <- ComputeGraphStatistics(toy.graph)
  Unodes <- c(1:3, 7:8)
  Vnodes <- setdiff(1:graph.stats$N, Unodes)
  p.score.obj <- .GetPscoreObjectBipartite(
      comm=7:10, graph.stats$adj, toy.graph, graph.stats$d,
      graph.stats$N, graph.stats$m, Unodes, Vnodes)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.01818182)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.01818182)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 7)
})

test_that("p-score function works on small imbalanced barbell graph.", {
  # Compute p-scores on a 2x2 bipartite clique next to a 3x3 bipartite clique
  # that has one cross-edge.
  # Let U = {1, 2, 3, 7, 8} and V = {4, 5, 6, 9, 10}
  # The cliques are {1, 2, 3, 4, 5, 6} and {7, 8, 9, 10}.
  # The extra edge connects 3 and 9.
  # Let UC = {7, 8} and VC = {9, 10}, UB = {1, 2, 3}, VB = {4, 5, 6}.
  # Considering U-side tests, some constant values:
  #   VB = {4, 5, 6} has 3 * 3 = 9 stubs.
  #   VC = {9, 10} has one outgoing stubs.
  #   Removing any node from UC = {7, 8} adds 2 outgoing stubs to VC.
  #   Each node in UC also has 2 outgoing stubs, and two observed white ball connections.
  #   So the p-scores should be identically:
  #     phyper(2 - 1, 3, 9, 2, lower.tail=FALSE) = 0.04545455
  #   ...with lower bounds of zero.
  # Considering V-side tests, some constant values:
  #   VB = {1, 2, 3} has 3 * 3 + 1 = 10 stubs.
  #   VC = {7, 8} has zero outgoing stubs.
  #   Removing any node from UC = {9, 10} adds 2 outgoing stubs to VC.
  #   Each node in UC has two observed white ball connections, but 9 has 3 outgoing stubs.
  #   So the p-scores should be identically:
  #     phyper(2 - 1, 2, 10, 3, lower.tail=FALSE) = 0.04545455
  #     phyper(2 - 1, 2, 10, 2, lower.tail=FALSE) = 0.01515152
  #   ...each having a lower bound of 0, even though one has degree 3, since in both
  #   cases there are only two white balls.
  # The worst node is a tie between 7, 8, 9.
  toy.graph <- AddClique(2, 2, AddClique(3, 3))
  toy.graph <- rbind(toy.graph, matrix(c(3, 9), ncol=2))
  graph.stats <- ComputeGraphStatistics(toy.graph)
  Unodes <- c(1:3, 7:8)
  Vnodes <- setdiff(1:graph.stats$N, Unodes)
  p.score.obj <- .GetPscoreObjectBipartite(
      comm=7:10, graph.stats$adj, toy.graph, graph.stats$d,
      graph.stats$N, graph.stats$m, Unodes, Vnodes)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.04545455)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.04545455)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 7)
})

test_that("FOCS correctly scores (3x3)-clique in imbalanced barbell graph.", {
  # Feed FOCS the imbalanced barbell graph, the community C = {7, 8, 9, 10}, and
  # p parameter 0.5, so that it will test the first two worst nodes.
  # Also, we set nrand=1 so we can test the randomized scoring.

  # Following the reasoning above, the worst and 2nd-worst nodes (7 & 8) have
  # p-score upper/lower bounds of identically 0.04545455/0.0. Thus the first round
  # of testing should perform the following calculation:
  #   p <- 0.5
  #   nrand <- 1
  #   k <- round(4 * p)
  #   fMat <- matrix(0, nrow=nrand, ncol=k)
  #   RUN_SEED <- 543212345
  #   set.seed(RUN_SEED)
  #
  #   # First round of testing.
  #   p.rand <- runif(2, c(0.0, 0.0), c(0.04545455, 0.04545455))
  #   N <- 10
  #   NC <- 4
  #   if (p.rand[1] > p.rand[2]) {
  #     test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  #   } else {
  #     test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  #   }
  #   fMat[1, 1] <- test.score

  # After 7 is removed, C = {8, 9, 10}, UC = {8}, VC = {9, 10},
  # UB = {1, 2, 3, 7}, VB = {4, 5, 6}. Considering U-side tests,
  # VB still has 9 stubs. VC has 3 outgoing stubs. Node 8 has zero outgoing stubs,
  # so removing it just adds 2 outgoing stubs to VC (now all its stubs are outgoing).
  # Thus node 8 has a p-score of
  #   phyper(2 - 1, 5, 9, 2, lower.tail=FALSE) = 0.1098901,
  # and a lower bound of zero.
  # Considering V-side tests, UB has 3 * 3 + 1 + 2 = 12 stubs. UC has no
  # outgoing stubs. Removing node 9 and 10 from VC each adds an outgoing stub
  # to UC. So the p-scores are
  #   phyper(1 - 1, 1, 12, 3, lower.tail=FALSE) = 0.2307692
  #   phyper(1 - 1, 1, 12, 2, lower.tail=FALSE) = 0.1538462
  # with lower bounds of zero since there is only one white ball.

  # Thus the second round of testing should perform the following calculation:
  #   p.rand <- runif(2, c(0.0, 0.0), c(0.2307692, 0.1538462))
  #   N <- 10
  #   NC <- 3
  #   if (p.rand[1] > p.rand[2]) {
  #     test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  #   } else {
  #     test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  #   }
  #   fMat[1, 2] <- test.score

  # Testing p-scores for second round of testing (first round in previous test case).
  toy.graph <- AddClique(2, 2, AddClique(3, 3))
  toy.graph <- rbind(toy.graph, matrix(c(3, 9), ncol=2))
  graph.stats <- ComputeGraphStatistics(toy.graph)
  Unodes <- c(1:3, 7:8)
  Vnodes <- setdiff(1:graph.stats$N, Unodes)
  p.score.obj <- .GetPscoreObjectBipartite(
      comm=8:10, graph.stats$adj, toy.graph, graph.stats$d,
      graph.stats$N, graph.stats$m, Unodes, Vnodes)

  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.2307692)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.1538462)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 9)

  # Testing rounds of p-score randomization.
  p <- 0.5
  nrand <- 1
  k <- round(4 * p)
  fMat <- matrix(0, nrow=nrand, ncol=k)
  RUN_SEED <- 543212345
  set.seed(RUN_SEED)

  # First round of testing.
  p.rand <- runif(2, c(0.0, 0.0), c(0.04545455, 0.04545455))
  N <- 10
  NC <- 4
  if (p.rand[1] > p.rand[2]) {
    test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  } else {
    test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  }
  fMat[1, 1] <- test.score

  p.rand <- runif(2, c(0.0, 0.0), c(0.2307692, 0.1538462))
  N <- 10
  NC <- 3
  if (p.rand[1] > p.rand[2]) {
    test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  } else {
    test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  }
  fMat[1, 2] <- test.score

  expect_equal(signif(median(as.vector(apply(fMat, 1, min))), 7),
               0.05405928)
  set.seed(RUN_SEED)
  expect_equal(.focs(7:10, graph.stats$adj, toy.graph, graph.stats$d,
                     graph.stats$N, graph.stats$m, p=0.25, nrand=1, Unodes),
               signif(median(as.vector(apply(fMat, 1, min))), 7))
})

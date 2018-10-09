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
 
context("FOCS Backend Check")
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

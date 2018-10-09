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

MakeToyGraph <- function () {
# This function makes a 5-node graph for testing FOCS:
#
#  e
#   \
#    \
#     2
#      \
#       \
#  a--1--b
#  |    /|
#  |   / |
#  1  2  3
#  | /   |
#  |/    |
#  c--3--d
  edges <- matrix(c(1, 2, 1,
                    1, 3, 1,
                    2, 3, 2,
                    2, 4, 3,
                    2, 5, 2,
                    3, 4, 3),
                    ncol=3, byrow=TRUE)
  return(edges)
}

AddClique <- function (k, edges=NULL) {

  if (is.null(edges)) {
    N <- 0
  } else {
    N <- max(edges)
  }

  new.edges <- matrix(0, ncol=2, nrow=k * (k - 1) / 2)
  counter <- 1
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
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

test_that("p-score function works on small toy graph.", {
  # Compute p-scores for community C = {b, c, d} in the toy graph.
  # .GetPscoreObject will compute the worst node in C, and the worst
  # and 2nd-worst p-scores (and their lower bounds). To check them,
  # here is how to manually compute the 3 p-scores for C:
  #  0. Here are some constant values:
  #       C' has 4 stubs in the model.
  #       C has 4 outgoing stubs in the model.
  #  1. For node b:
  #       b has 5 edges into the community, and 3 leaving it.
  #       So, if b is removed, C now has 4 + 5 - 3 = 6 outgoing stubs.
  #       Also, b has 8 stubs overall.
  #       Therefore, the pscore of b is:
  #         phyper(5 - 1, 6, 4, 8, lower.tail=FALSE) = 0.6666667 
  #       
  #  2. Similarly, for node c:
  #         phyper(5 - 1, 8, 4, 6, lower.tail=FALSE) = 0.2727273
  #
  #  3. Similarly, for node d:
  #         phyper(6 - 1, 10, 4, 6, lower.tail=FALSE) = 0.06993007 
  # Thus, the worst node is b, and the 2nd worst is c. The lower bounds
  # of those p-scores are:
  #         phyper(5, 6, 4, 8, lower.tail=FALSE) = 0.1333333 (for b)
  #         phyper(5, 8, 4, 6, lower.tail=FALSE) = 0.03030303 (for c)
  # The first p-score is largest, so worst = 2.
  toy.graph <- MakeToyGraph() 
  graph.stats <- ComputeGraphStatistics(toy.graph)
  p.score.obj <- .GetPscoreObject(comm=2:4, graph.stats$adj, toy.graph,
                                  graph.stats$d, graph.stats$N, graph.stats$m)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.6666667)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.2727273)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0.1333333)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0.03030303)
  expect_equal(p.score.obj$worst, 2)
})

test_that("p-score function works on small clique graph.", {
  # Compute p-scores on a 3-clique next to a 5-clique.
  # Some constant values:
  #   C' has 5 * 4 = 20 stubs.
  #   C  has 0 outgoing stubs.
  #   Removing any node adds 2 outgoing stubs to C.
  #   Every node also has 2 outgoing stubs.
  #   So the p-scores should be identically:
  #     phyper(2 - 1, 2, 20, 2, lower.tail=FALSE) = 0.004329004
  #   ...with lower bounds of zero.
  # The p-scores are tied, so worst = 6. 
  toy.graph <- AddClique(3, AddClique(5)) 
  graph.stats <- ComputeGraphStatistics(toy.graph)
  p.score.obj <- .GetPscoreObject(comm=6:8, graph.stats$adj, toy.graph,
                                  graph.stats$d, graph.stats$N, graph.stats$m)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.004329004)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.004329004)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 6)
})

test_that("p-score function works on small imbalanced barbell graph.", {
  # Compute p-scores on a 3-clique next to a 5-clique, connected by an edge.
  # Some constant values:
  #   C' has 5 * 4 + 1 = 21 stubs.
  #   C  has 1 outgoing stubs.
  #   Removing node 6 from {6, 7, 8} adds 2 outgoing stubs to C but removes 1.
  #   Removing either 7 or 8 adds 2 outgoing stubs to C.
  #   Every node also has 2 outgoing stubs.
  #   So the p-score for 6 should be:
  #     phyper(2 - 1, 2, 21, 3, lower.tail=FALSE) = 0.01185771
  #   And the p-score for 7 and 8 should be:
  #     phyper(2 - 1, 3, 21, 2, lower.tail=FALSE) = 0.01086957
  #   ...with lower bounds of zero for both, since:
  #     1. for node 6 you can only get 2 white balls (even though it has degree 3)
  #     2. for node 7/8 you can get 3 white balls but each only has degree 2.
  # The p-scores are tied, so worst = 6.
  toy.graph <- AddClique(3, AddClique(5))
  toy.graph <- rbind(toy.graph, matrix(c(5, 6), ncol=2))
  graph.stats <- ComputeGraphStatistics(toy.graph)
  p.score.obj <- .GetPscoreObject(comm=6:8, graph.stats$adj, toy.graph,
                                  graph.stats$d, graph.stats$N, graph.stats$m)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.01185771)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.01086957)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 6)
})

test_that("FOCS correctly scores 5-clique in imbalanced barbell graph.", {

  # Feed FOCS the imbalanced barbell graph, the community {1, 2, 3, 4, 5}, and
  # p parameter 0.5, so that it will test the first two worst nodes. 
  # Also, we set nrand=1 so we can test the randomized scoring.

  # At first we have C' with 3 * 2 + 1 = 7 stubs, and C with 1 outgoing stub.
  # Removing node 5 from C adds 4 outgoing stubs to C but removes 1.
  # Removing any other nodes in C adds 4 outgoing stubs to C.
  # So the p-score for 5 should be:
  #   phyper(4 - 1, 4, 7, 5, lower.tail=FALSE) = 0.01515152.
  # and for {1...4}:
  #   phyper(4 - 1, 5, 7, 4, lower.tail=FALSE) = 0.01010101.
  # ...with lower bounds of zero for both, due to reasoning in previous test.
  # So the first set of worst nodes are 5 and 1, and 5 will be removed.

  # If the seed is set at 543212345, the .fScore function should perform a
  # calculation identical to these commands:
  #   p <- 0.5
  #   nrand <- 1
  #   k <- round(5 * p)
  #   fMat <- matrix(0, nrow=nrand, ncol=k)
  #   RUN_SEED <- 543212345
  #   set.seed(RUN_SEED)
  #
  #   # First round of testing.
  #   p.rand <- runif(2, c(0.0, 0.0), c(0.01515152, 0.01010101))
  #   N <- 8
  #   NC <- 5
  #   if (p.rand[1] > p.rand[2]) {
  #     test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)  
  #   } else {
  #     test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  #   }
  #   fMat[1, 1] <- test.score
  #
  # After 5 is removed, C = {1, 2, 3, 4}, C' nodes have the following stub counts:
  #   5: 5, 6: 3, 7: 2, 7: 2,
  # so 12 stubs total. C has 4 outgoing stubs to node 5 in C'. Removing any node in C
  # adds 3 outgoing stubs but removes 1 (the one it had to node 5). So p-scores are:
  #   phyper(3 - 1, 6, 12, 4, lower.tail=FALSE) = 0.08333333. 
  # with lower bounds:
  #   phyper(3, 6, 12, 4, lower.tail=FALSE) = 0.004901961.
  #
  # With the worst node removed, .focs should call .fScore again (i = 2) with the above
  # upper/lower bounds:
  #
  #   # Second round of testing.
  #   p.rand <- runif(2, c(0.0, 0.0), c(0.08333333, 0.004901961))
  #   N <- 8
  #   NC <- 4
  #   if (p.rand[1] > p.rand[2]) {
  #     test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  #   } else {
  #     test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  #   }
  #   fMat[1, 2] <- test.score
  #
  toy.graph <- AddClique(3, AddClique(5))
  toy.graph <- rbind(toy.graph, matrix(c(5, 6), ncol=2))
  graph.stats <- ComputeGraphStatistics(toy.graph)

  p.score.obj <- .GetPscoreObject(comm=1:5, graph.stats$adj, toy.graph,
                                  graph.stats$d, graph.stats$N, graph.stats$m)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.01515152)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.01010101)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0)
  expect_equal(p.score.obj$worst, 5)
  
  p.score.obj <- .GetPscoreObject(comm=1:4, graph.stats$adj, toy.graph,
                                  graph.stats$d, graph.stats$N, graph.stats$m)
  expect_equal(signif(p.score.obj$ps.upper[1], 7), 0.08333333)
  expect_equal(signif(p.score.obj$ps.lower[1], 7), 0.004901961)
  expect_equal(signif(p.score.obj$ps.upper[2], 7), 0.08333333)
  expect_equal(signif(p.score.obj$ps.lower[2], 7), 0.004901961)
  expect_equal(p.score.obj$worst, 1) 

  p <- 0.5
  nrand <- 1
  k <- round(5 * p)
  fMat <- matrix(0, nrow=nrand, ncol=k)
  RUN_SEED <- 543212345
  set.seed(RUN_SEED)
  # First round of testing: note this does not merely duplicate code,
  # we are (1) hard-coding the computed p-scores, and 
  # (2) manually executing the two loop rounds.
  p.rand <- runif(2, c(0.0, 0.0), c(0.01515152, 0.01010101))
  N <- 8
  NC <- 5
  if (p.rand[1] > p.rand[2]) {
    test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  } else {
    test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  }
  fMat[1, 1] <- test.score

  # Second round of testing: note this does not merely duplicate code,
  # we are (1) hard-coding the computed p-scores, and 
  # (2) manually executing the two loop rounds.
  p.rand <- runif(2, c(0.0, 0.0), c(0.08333333, 0.004901961))
  N <- 8
  NC <- 4
  if (p.rand[1] > p.rand[2]) {
    test.score <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
  } else {
    test.score <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
  }
  fMat[1, 2] <- test.score
  
  expect_equal(signif(median(as.vector(apply(fMat, 1, min))), 7),
               0.004881726)
  set.seed(RUN_SEED)
  expect_equal(.focs(1:5, graph.stats$adj, toy.graph, graph.stats$d, 
                     graph.stats$N, graph.stats$m, p=0.25, nrand=1),
               signif(median(as.vector(apply(fMat, 1, min))), 7))
})

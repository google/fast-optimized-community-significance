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
# 
#' Compute significance of a single community.
#'
#' @param comm (integer) The community.
#' @param adjG (dgCMatrix) The adjacency matrix of the network.
#' @param edges (integer matrix) The edgelist of the network.
#' @param d (integer) The degree vector of the network.
#' @param N (integer) The number of nodes in the network.
#' @param m (integer) Twice the number of edges in the network.
#' @param nrand (integer, default: 100) Number of random choices of the p-scores.
.fScore <- function (comm, adjG, edges, d, N, m, nrand) { 

  # Compute basic values.
  d_u_Cs <- Matrix::colSums(adjG[comm, ])
  d_u_Cprimes <- d - d_u_Cs
  d_C <- sum(d[comm])
  d_C_C <- sum(d_u_Cs[comm])
  d_C_Cprime <- d_C - d_C_C
  d_Cprime_Cprime <- sum(d) - d_C - d_C_Cprime
  NC <- length(comm)
  
  # Compute p-scores.
  pScores <- phyper(q=d_u_Cs[comm] - 1,
                    m=d_C_Cprime - d_u_Cprimes[comm],  
                      # white balls: edges coming out of C.
                      # subtraction removes out-degree of community nodes.
                    n=d_Cprime_Cprime + d_u_Cprimes[comm], 
                      # black balls: edges totally outside C.
                      # addition adds in out-degree of community nodes.
                    k=d[comm],
                    lower.tail=FALSE)

  # Find worst node indexs.
  pScoresOrd <- order(pScores, decreasing=TRUE)
  w <- comm[pScoresOrd[1]]
  v <- comm[pScoresOrd[2]]


  # Compute lower/upper pScore bounds for w and v with w removed.
  pU <- phyper(d_u_Cs[c(w, v)] - 1, 
               d_C_Cprime - d_u_Cprimes[w],
               d_Cprime_Cprime + d_u_Cprimes[w],
               d[c(w, v)], lower.tail=FALSE)
  pL <- phyper(d_u_Cs[c(w, v)], 
               d_C_Cprime - d_u_Cprimes[w],
               d_Cprime_Cprime + d_u_Cprimes[w],
               d[c(w, v)], lower.tail=FALSE)

  # Get random pscores and compute test.
  cs <- numeric(nrand)
  for (counter in 1:nrand) {    
    pRand <- runif(2, pL, pU)
    if (pRand[1] > pRand[2]) {
      cs[counter] <- 1 - ((1 - pRand[1]) / (1 - pRand[2]))^(N - NC + 1)
    } else {
      cs[counter] <- 1 - ((1 - pRand[2]) / (1 - pRand[1]))^(N - NC + 1)
    }
  }

  # Return.
  return(list(scores=cs, worst=w))
}

#' FOCS backend.
#'
#' @param comm (integer) The community.
#' @param adjG (dgCMatrix) The adjacency matrix of the network.
#' @param edges (integer matrix) The edgelist of the network.
#' @param d (integer) The degree vector of the network.
#' @param N (integer) The number of nodes in the network.
#' @param m (integer) Twice the number of edges in the network.
#' @param p (double, default: 0.25) Proportion of community to use as the border.
#' @param nrand (integer, default: 100) Number of random choices of the p-scores.
.focs <- function (comm, adjG, edges, d, N, m, p=0.25, nrand=100) { 
  if (length(comm) < 3) {
    return(1)
  }
  
  nodes <- 1:N
  k <- counter <- ceiling(length(comm) * p)
  fMat <- matrix(0, nrow=nrand, ncol=k)

  for (i in 1:k) {
    fScoreList <- .fScore(comm, adjG, edges, d, N, m, nrand)
    fMat[, i] <- fScoreList$scores
    comm <- setdiff(comm, fScoreList$worst)
  }
  minFScores <- as.vector(apply(fMat, 1, min))
  return(median(minFScores))
}

#' Run FOCS significance computation on a list of communities.
#'
#' @param community_list (list of int vectors) List of communities to compute significance of.
#' @param edges (2/3-column integer matrix) Rows are edge indexes and (optionally) edge weight.
#' @param p (double, default: 0.25) Proportion of community to use as the border.
#' @param nrand (integer, default: 100) Number of random choices of the p-scores.
#' @export
FOCS <- function (community_list, edges, p=0.25, nrand=100) {
  N <- max(edges[, c(1:2)])
  if (ncol(edges) == 2) {
    values <- rep(1, nrow(edges))
  } else {
    values <- edges[, 3]
  }
  adj <- Matrix::sparseMatrix(i=edges[, 1], j=edges[, 2],
                              x=values, dims=c(N, N), symmetric=TRUE)
  d <- Matrix::colSums(adj)
  N <- length(d)
  m <- nrow(edges)
  scores <- lapply(community_list, .focs, adj, edges, d, N, m, p, nrand)
  return(scores)
}

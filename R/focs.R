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

#' Get Pscores for a community.
#'
#' @param observations (integer vector) Number of observed white balls.
#' @param cross.edges (integer) number of possible white balls.
#' @param outside.edges (integer) number of possible black balls.
#' @param draws (integer vector of length(observations)) Possible draws.
#' @param cross.edge.correction (integer or integer vector of length(observations)) Additions to #whites after removing nodes.
.GetPscores <- function (observations, cross.edges, outside.edges, draws,
                         cross.edge.correction=0) {
  return(phyper(q=observations - 1,
                m=cross.edges + cross.edge.correction,
                # white balls: edges coming out of C, including those toward 
                #              removed nodes, if applicable.
                n=outside.edges,
                # black balls: edges within or coming out of C'.
                k=draws,
                lower.tail=FALSE))
}

#' Get pscores for worst & 2nd-worst node, and index of worst node.
#
#' @param comm (integer) The community.
#' @param adjG (dgCMatrix) The adjacency matrix of the network.
#' @param edges (integer matrix) The edgelist of the network.
#' @param d (integer) The degree vector of the network.
#' @param N (integer) The number of nodes in the network.
#' @param m (integer) Twice the number of edges in the network.

.GetPscoreObject <- function (comm, adjG, edges, d, N, m) {

  # Compute basic values.
  d_u_Cs <- Matrix::colSums(adjG[comm, ])
  d_u_Cprimes <- d - d_u_Cs
  d_C <- sum(d[comm])
  d_Cprime <- sum(d) - d_C
  d_C_C <- sum(d_u_Cs[comm])
  d_C_Cprime <- d_C - d_C_C
  d_Cprime_Cprime <- sum(d) - d_C - d_C_Cprime
  NC <- length(comm)
  
  # Compute p-scores.
  p.scores <- .GetPscores(observations=d_u_Cs[comm], cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[comm],
                          cross.edge.correction=d_u_Cs[comm] - d_u_Cprimes[comm])

  # Find worst node indexs.
  p.scores.ord <- order(p.scores, decreasing=TRUE)
  w <- comm[p.scores.ord[1]]
  v <- comm[p.scores.ord[2]]

  # Re-compute p-scores for w & v, and their upper bounds.
  ps.upper <- .GetPscores(observations=d_u_Cs[c(w, v)], cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[c(w, v)],
                          cross.edge.correction=d_u_Cs[c(w, v)] - d_u_Cprimes[c(w, v)])
  ps.lower <- .GetPscores(observations=d_u_Cs[c(w, v)] + 1, cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[c(w, v)],
                          cross.edge.correction=d_u_Cs[c(w, v)] - d_u_Cprimes[c(w, v)])

  return(list(worst=w, ps.upper=ps.upper, ps.lower=ps.lower)) 
}

#' Compute significance of a single community.
#'
#' params inherit from .GetPscoreObject.
.fScore <- function (comm, adjG, edges, d, N, m, nrand) { 
 
 ps.object <- .GetPscoreObject(comm, adjG, edges, d, N, m, nrand)
 NC <- length(comm)

 # Get random pscores and compute test.
  cs <- numeric(nrand)
  for (counter in 1:nrand) {    
    p.rand <- runif(2, ps.object$ps.lower, ps.object$ps.upper)
    if (p.rand[1] > p.rand[2]) {
      cs[counter] <- 1 - ((1 - p.rand[1]) / (1 - p.rand[2]))^(N - NC + 1)
    } else {
      cs[counter] <- 1 - ((1 - p.rand[2]) / (1 - p.rand[1]))^(N - NC + 1)
    }
  }

  # Return.
  return(list(scores=cs, worst=ps.object$worst))
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
  k <- round(length(comm) * p)
  if (k == 0) {
    k <- 1
  }
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
  m <- sum(d) / 2
  scores <- lapply(community_list, .focs, adj, edges, d, N, m, p, nrand)
  return(scores)
}

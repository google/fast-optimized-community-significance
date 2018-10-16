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
  d_u_Cs <- Matrix::rowSums(adjG[comm, comm])
  d_u_Cprimes <- d[comm] - d_u_Cs
  d_C <- sum(d[comm])
  d_Cprime <- sum(d) - d_C
  d_C_C <- sum(d_u_Cs)
  d_C_Cprime <- d_C - d_C_C

  # Compute p-scores.
  p.scores <- .GetPscores(observations=d_u_Cs, cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[comm],
                          cross.edge.correction=d_u_Cs - d_u_Cprimes)

  # Find worst node indexs.
  p.scores.ord <- order(p.scores, decreasing=TRUE)
  w <- p.scores.ord[1]
  v <- p.scores.ord[2]

  # Re-compute p-scores for w & v, and their upper bounds.
  ps.upper <- .GetPscores(observations=d_u_Cs[c(w, v)], cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[comm[c(w, v)]],
                          cross.edge.correction=d_u_Cs[c(w, v)] - d_u_Cprimes[c(w, v)])
  ps.lower <- .GetPscores(observations=d_u_Cs[c(w, v)] + 1, cross.edges=d_C_Cprime,
                          outside.edges=d_Cprime, draws=d[comm[c(w, v)]],
                          cross.edge.correction=d_u_Cs[c(w, v)] - d_u_Cprimes[c(w, v)])

  return(list(worst=comm[w], ps.upper=ps.upper, ps.lower=ps.lower))
}

#' Get pscores for worst & 2nd-worst node, and index of worst node (bipartite version).
#
#' @param Unodes (integer vector) nodes on the U side.
#' @param Vnodes (integer vector) nodes on the V side.
#' other params inherit from .GetPscoreObject.
.GetPscoreObjectBipartite <- function (comm, adjG, edges, d, N, m, Unodes, Vnodes) {
  if (mean(comm %in% Unodes) %in% c(0, 1)) {
    return(list(worst=sample(comm, 1), ps.upper=rep(1, 2), ps.lower=rep(1, 2)))
  }

  # Compute basic values.
  d_u_Cs <- Matrix::rowSums(adjG[comm, comm])
  d_u_Cprimes <- d[comm] - d_u_Cs

  # Compute p.scores for each side.
  sides <- list(Unodes, Vnodes)
  p.scores <- numeric(length(comm))
  for (i in 1:2) {
    comm.side <- intersect(sides[[i]], comm)
    comm.aside <- intersect(sides[[3 - i]], comm)
    comm.side.match <- match(comm.side, comm)
    d_Cprime <- m - sum(d[comm.aside])
    d_C_Cprime <- sum(d[comm.aside]) - sum(d_u_Cs[comm.side.match])
    p.scores[comm.side.match] <- .GetPscores(
        observations=d_u_Cs[comm.side.match], cross.edges=d_C_Cprime,
        outside.edges=d_Cprime, draws=d[comm.side],
        cross.edge.correction=d_u_Cs[comm.side.match])
  }

  # Find worst node indexs.
  p.scores.ord <- order(p.scores, decreasing=TRUE)
  w <- p.scores.ord[1]
  v <- p.scores.ord[2]

  # Re-compute p-scores for w & v, and their upper bounds.
  ps.upper <- ps.lower <- numeric(2)
  choice.nodes <- c(w, v)
  for (i in 1:2) {
    for (j in 1:2) {
      comm.side <- intersect(sides[[i]], comm)
      comm.aside <- intersect(sides[[3 - i]], comm)
      comm.side.match <- match(comm.side, comm)
      if (choice.nodes[j] %in% comm.side.match) {
        d_Cprime <- m - sum(d[comm.aside])
        d_C_Cprime <- sum(d[comm.aside]) - sum(d_u_Cs[comm.side.match])
        ps.upper[j] <- .GetPscores(observations=d_u_Cs[choice.nodes[j]], cross.edges=d_C_Cprime,
                                   outside.edges=d_Cprime, draws=d[comm[choice.nodes[j]]],
                                   cross.edge.correction=d_u_Cs[choice.nodes[j]])
        ps.lower[j] <- .GetPscores(observations=d_u_Cs[choice.nodes[j]] + 1, cross.edges=d_C_Cprime,
                                   outside.edges=d_Cprime, draws=d[comm[choice.nodes[j]]],
                                   cross.edge.correction=d_u_Cs[choice.nodes[j]])
      }
    }
  }
  return(list(worst=comm[w], ps.upper=ps.upper, ps.lower=ps.lower))
}

#' Helper function to compute score to high precision.
score.computer <- function(p1, p2, N, NC) {
  n <- (N - NC + 1)
  if (p1 * n < 0.01 && p2 * n < 0.01) {
    num <- 1 - n * p1 + n * (n - 1) * p1^2 / 2
    den <- 1 - n * p2 + n * (n - 1) * p2^2 / 2
    return(1 - num / den)
  } else {
    return(1 - ((1 - p1) / (1 - p2))^n)
  }
}

#' Compute significance of a single community.
#'
#' params inherit from .GetPscoreObject.
.fScore <- function (comm, adjG, edges, d, N, m, nrand,
                     Unodes=NULL, Vnodes=NULL) {

  if (is.null(Unodes)) {
    ps.object <- .GetPscoreObject(comm, adjG, edges, d, N, m)
  } else {
    ps.object <- .GetPscoreObjectBipartite(comm, adjG, edges, d, N, m,
                                           Unodes, Vnodes)
  }
  NC <- length(comm)

  # Get random pscores and compute test.
  cs <- numeric(nrand)
  for (counter in 1:nrand) {
    p.rand <- runif(2, ps.object$ps.lower, ps.object$ps.upper)
    if (p.rand[1] > p.rand[2]) {
      cs[counter] <- score.computer(p.rand[1], p.rand[2], N, NC)
    } else {
      cs[counter] <- score.computer(p.rand[2], p.rand[1], N, NC)
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
#' @param Unodes (integer vector, default: NULL) Node list of U-side nodes. If supplied, bipartite analysis is assumed.
.focs <- function (comm, adjG, edges, d, N, m, p=0.25, nrand=100, Unodes=NULL) {
  if (length(comm) < 3) {
    return(1)
  }

  if (is.null(Unodes)) {
    Vnodes <- NULL
  } else {
    Vnodes <- setdiff(1:N, Unodes)
  }
  k <- round(length(comm) * p)
  if (k == 0) {
    k <- 1
  }
  fMat <- matrix(0, nrow=nrand, ncol=k)

  for (i in 1:k) {
    fScoreList <- .fScore(comm, adjG, edges, d, N, m, nrand, Unodes, Vnodes)
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
#' @param Unodes (integer vector, default: NULL) Node list of U-side nodes. If supplied, bipartite analysis is assumed.
#' @export
FOCS <- function (community_list, edges, p=0.25, nrand=100, Unodes=NULL) {
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
  scores <- lapply(community_list, .focs, adj, edges, d, N, m, p, nrand, Unodes)
  return(scores)
}


#' Compute statistical significance of a single community.
#'
#' @param comm (integer) The community.
#' @param G (igraph) The network in igraph format.
#' @param adjG (dgCMatrix) The adjacency matrix of the network.
#' @param edges (integer matrix) The edgelist of the network.
#' @param d (integer) The degree vector of the network.
#' @param N (integer) The number of nodes in the network.
#' @param m (integer) Twice the number of edges in the network.
#' @param version (integer, default: 3) Version of the method to use.
#' @param borderp (double, default: 0.25) Proportion of community to use as the border.
#' @param nmcsims (integer, default: 100) Number of monte carlo simulations of p-values to produce.
#' @param return_mean (bool, default: FALSE) Whether to return mean or median of mc p-values.
.SOCS <- function (comm, adjG, edges, d, N, m, borderp=0.25, nmcsims=100) { 

  if (length(comm) < 3)
    return(1)
  
  # Computing basic values
  nodes <- 1:N
  commInDegrees <- Matrix::colSums(adjG[comm, ])
  commOutDegrees <- d - commInDegrees
  commInDegT <- sum(commInDegrees[comm])
  commDegT <- sum(d[comm])
  commOutDegT <- commDegT - commInDegT
  extDeg <- 2 * m - commDegT
  
  # Getting r-scores
  commOutDegT2 <- commOutDegT - commOutDegrees[comm] + commInDegrees[comm]
  commDegT2 <- commDegT - d[comm]
  extDeg2 <- 2 * m - commDegT2
  pScores <- phyper(commInDegrees[comm] - 1, commOutDegT2, 
                     extDeg2 - commOutDegT2, d[comm], lower.tail=FALSE)
  worst_node <- comm[which.max(pScores)]

  cs <- numeric(nmcsims)
  commc <- setdiff(nodes, comm)
  pScoresU <- phyper(commInDegrees[comm] - 1, commOutDegT2, 
                      extDeg2 - commOutDegT2, d[comm], lower.tail=FALSE)
  pScoresL <- phyper(commInDegrees[comm], commOutDegT2,
                      extDeg2 - commOutDegT2, d[comm], lower.tail=FALSE)
    
  for (counter in 1:nmcsims) {    
    # Getting random SOCS scores
    pScoresRand <- runif(length(comm), pScoresL, pScoresU)
    k <- ceiling(length(comm) * borderp)
    scs <- numeric(k)
    pScoresSort <- sort(pScoresRand, decreasing=TRUE)
    invP <- (1 - pScoresSort[1:k]) / (1 - pScoresSort[2:(k + 1)])
    scs <- 1 - invP^(N - length(comm) + 2:(k + 1))
    cs[counter] <- min(scs) * k      
  }

  return(median(cs))      
}

#' Compute statistical significance of a single community.
#'
#' @param community_list (list of integers) List of communities to compute significance of.
#' @param network (igraph) The network in igraph format.
#' @param borderp (double, default: 0.25) Proportion of community to use as the border.
#' @param nmcsims (integer, default: 100) Number of monte carlo simulations of p-values to produce.
#' @export
SOCS <- function (community_list, network, borderp=0.25, nmcsims=100) {

  # Computing network statistics
  adj <- get.adjacency(network)
  d <- degree(network)
  N <- length(V(network))
  edges <- get.edgelist(network)
  m <- nrow(edges)

  K <- length(community_list)
  pvals <- numeric(K)
  for (i in 1:K) {
    pvals[i] <- MtSigCommSingle(community_list[[i]], adj, edges, d, N, m,
                                borderp=borderp, nmcsims=nmcscims)
  }

  return(pvals)
}

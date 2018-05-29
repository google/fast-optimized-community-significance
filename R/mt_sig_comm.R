#' Compute statistical significance of a single community.
#'
#' @param group_s (integer) The community.
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
MtSigCommSingle <- function (group_s, G, adjG, edges, d, N, m, version=3, 
                             borderp=0.25, nmcsims=100, return_mean=FALSE) { 

  if (length(group_s) < 3)
    return(1)
  
  # Computing basic values
  nodes <- 1:N
  kints <- Matrix::colSums(adjG[group_s, ])
  kexts <- d - kints
  mC_int <- sum(kints[group_s])
  mC <- sum(d[group_s])
  mC_ext <- mC - mC_int
  mstar <- m - mC
  
  # Getting r-scores
  mC_ext_new <- mC_ext - kexts[group_s] + kints[group_s]
  mC_new <- mC - d[group_s]
  mstar_new <- m - mC_new
  r_scores <- phyper(kints[group_s] - 1, mC_ext_new, mstar_new - mC_ext_new,
                     d[group_s], lower.tail = FALSE)
  worst_node <- group_s[which.max(r_scores)]
  kw <- kints[worst_node]
  degree_ref <- d[worst_node]
  
  if (version > 1) {
    
    cs <- numeric(nmcsims)
    group_sc <- setdiff(nodes, group_s)
    r_scoresU <- phyper(kints[group_s] - 1, mC_ext_new, mstar_new - mC_ext_new,
                        d[group_s], lower.tail = FALSE)
    r_scoresL <- phyper(kints[group_s],     mC_ext_new, mstar_new - mC_ext_new,
                        d[group_s], lower.tail = FALSE)
    
    for (counter in 1:nmcsims) {
      
      # Getting random r-scores
      rand_rscores <- runif(length(group_s), r_scoresL, r_scoresU)
      
      if (version == 2) {
        r1 <- sort(rand_rscores, decreasing = TRUE)[1]
        r2 <- sort(rand_rscores, decreasing = TRUE)[2]
        cs[counter] <- 1 - (1 - (r1 - r2) / (1 - r2))^(N - length(group_s) + 1)
      }
      
      if (version == 3) {
        
        ntest <- ceiling(length(group_s) * borderp)
        scs <- numeric(ntest)
        rscore_sort <- sort(rand_rscores, decreasing=TRUE)
        inv.p <- (1 - rscore_sort[1:ntest]) / (1 - rscore_sort[2:(ntest + 1)])
        scs <- 1 - inv.p^(N - length(group_s) + 2:(ntest + 1))
        cs[counter] <- min(scs) * ntest
        
      }
    }
    
    
    if (return_mean) {
      return(mean(cs))
    } else {
      return(median(cs))
    }
  }
  
  
  if (FALSE) {
    
    # Doing weird step that came with the cb-signi package...so far has not
    # thrown an error, but it does on growing cliq graph
    r_threshold1 <- sort(r_scores, decreasing = TRUE)[2]
    worst_node2 <- group_s[order(r_scores, decreasing = TRUE)[2]]
    k2 <- kints[worst_node2]
    s2 <- d[worst_node2]
    r_threshold2 <- phyper(k2 - 1, mC_ext + 2 * k2 - degree_ref, 
                           mstar + s2 + degree_ref + 2 * k2 - mC_ext,
                           s2, lower.tail = FALSE)
    r_threshold <- min(r_threshold1, r_threshold2)
    kv_eq <- qhyper(r_threshold, mC_ext, mstar - mC_ext, 
                    degree_ref, lower.tail = FALSE)
    
  } else {
    
    # A more kosher-ized version of the weird step. Following paper exactly.
    # Doesn't use two different r thresholds.
    r_threshold <- sort(r_scores, decreasing = TRUE)[2]
    
    # Computing dist'n corresponding to worst node.
    mC_ext_now <- mC_ext - kexts[worst_node] + kints[worst_node]
    mC_now <- mC - d[worst_node]
    mstar_now <- m - mC_now
    
    # Some r_thresholds can be very small so we have to expand the whole dist'n
    all_qs <- phyper(0:degree_ref, mC_ext_now, mstar_now - mC_ext_now,
                     degree_ref, log.p = TRUE, lower.tail = FALSE)
    all_qs <- c(0, all_qs)
    #kv_eq <- max(which(all_qs >= log(r_threshold))) - 1
    
    # Correction for loss of precision
    kv_eq <- max(which(all_qs - log(r_threshold) >= -1e-10))
    
  }
  
  if (kv_eq < kw) {
    message("ERROR; report this message to johnpalowitch@gmail.com, thank you!")
    return(-1)
  }
  
  # Removing the worst node
  group_s <- setdiff(group_s, worst_node)
  
  # Computing mcout_averaged ---------------------------------------------------
  
  # Re-computing set values
  #tstrength <- sum(degree(G)) / 2
  kints <- colSums(adjG[group_s, ])
  kexts <- d - kints
  mC_int <- sum(kints[group_s])
  mC <- sum(d[group_s])
  mC_ext <- mC - mC_int
  mstar <- m - mC

  # Initializing loop
  mcout <- 0
  mcs <- numeric(0)
  avedevmax <- 0.1
  counter <- 0
  group_sc <- setdiff(nodes, group_s)
  
  while (TRUE) {
    counter <- counter + 1
    a1s <- phyper(kints[group_sc] - 1, mC_ext, mstar - mC_ext, d[group_sc])
    a2s <- phyper(kints[group_sc],     mC_ext, mstar - mC_ext, d[group_sc])
    rs <- runif(length(a1s), a1s, a2s)
    # Qs <- qhyper(rs, mC_ext, mstar - mC_ext, degree_ref)
    # Note: mC_ext, mstar - mC_ext, & degree_ref make dist'n of all_qs (L:133)
    Qs <- sapply(rs, function (r) max(which(all_qs - log(r) >= -1e-10)) - 1)
    mcs[counter] <- sum(Qs)
    avedev <- ifelse(length(mcs) > 10, sqrt(var(mcs) / length(mcs)), 1)
    if (avedev < avedevmax || counter > 500)
      break
  }
  
  
  mcout <- round(mean(mcs))
  #cat("--mcout is", mcout, "\n")
  cscore <- 0
	
	Ncc <- N - length(group_s)
  
	if (kw == 0) {
		cscore <- 1
	} else if (kv_eq * Ncc < mcout) {
    return(1)
	} else {
				
		#prob <- phyper(kw - 1, mcout, kv_eq * Ncc - mcout, kv_eq)
	  #cscore <- 1 - prob^Ncc
		qrob <- phyper(kw - 1, mcout, kv_eq * Ncc - mcout, kv_eq, lower.tail = FALSE)
		#cscore <- (1 - qrob)^Ncc
		#cscore <- 1 - (1 - qrob)^Ncc
		cscore <- qrob
		
	
	}
	
	#newrsU <- phyper(Qs - 1, mcout, kv_eq * Ncc - mcout, kv_eq, lower.tail = FALSE)
	#newrsL <- phyper(Qs, mcout, kv_eq * Ncc - mcout, kv_eq, lower.tail = FALSE)
	#newrs_rand <- runif(length(Qs), newrsL, newrsU)
	
	
	
	if(cscore<0)
		cscore=0
	
	
	return(cscore)
  
}

#' Compute statistical significance of a single community.
#'
#' @param community_list (list of integers) List of communities to compute significance of.
#' @param network (igraph) The network in igraph format.
#' @param borderp (double, default: 0.25) Proportion of community to use as the border.
#' @param nmcsims (integer, default: 100) Number of monte carlo simulations of p-values to produce.
#' @export
MtSigComm <- function (community_list, network, borderp=0.25, nmcsims=100) {

  # Computing network statistics
  adj <- get.adjacency(network)
  d <- degree(network)
  N <- length(V(network))
  m <- 2 * nrow(edges)
  edges <- get.edgelist(network)

  K <- length(community_list)
  pvals <- numeric(K)
  for (i in 1:K) {
    pvals[i] <- MtSigCommSingle(community_list[[i]], network, adj, edges, d, N, m,
                                borderp=borderp, nmcsims=nmcscims)
  }

  return(pvals)
}

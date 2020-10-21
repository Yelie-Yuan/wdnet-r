#' @importFrom stats weighted.mean
NULL

#' Directed assortativity coefficient
#' 
#' Compute the assortativity coefficient of a weighted and directed network.
#'
#' @usage
#' dw_assort(adj, type = c("out-in", "in-in", "out-out", "in-out"))
#'
#' @param adj is an adjacency matrix of an weighted and directed network.
#' @param type which type of assortativity coefficient to compute: "out-in" (default), 
#' "in-in", "out-out" or "in-out"?
#'
#' @return a scalar of assortativity coefficient
#'
#' @references
#' \itemize{
#' \item Foster, J.G., Foster, D.V., Grassberger, P. and Paczuski, M. (2010). Edge direction 
#' and the structure of networks. \emph{Proceedings of the National Academy of Sciences of the
#' United States}, 107(24), 10815--10820.
#' \item Yuan, Y. Zhang, P. and Yan, J. (2020+). Assortativity coefficients for 
#' weighted and directed networks
#' }
#'
#' @note 
#' When the adjacency matrix is binary (i.e., directed but unweighted networks), \code{dw_assort}
#' returns the assortativity coefficient proposed in Foster et al. (2010).
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' ## and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x) x*sample(3,1))
#' adj_ER <- matrix(weight_ER,20,20)
#' system.time(myassort <- dw_assort(adj_ER, type = "out-in"))
#' myassort
#' 
#' @export

dw_assort <- function(adj, type = c("out-in", "in-in", "out-out", "in-out")) {
  stopifnot(dim(adj)[1] == dim(adj)[2])
  ## determine the location of edges in the network
  in_str <- colSums(adj)
  out_str <- rowSums(adj)
  vert_from <- unlist(apply(adj, 2, function(x){which(x > 0)}))
  number_to <- apply(adj, 2, function(x){length(which(x > 0) == TRUE)})
  temp_to <- cbind(seq(1:dim(adj)[1]),number_to)
  vert_to <- rep(temp_to[,1],temp_to[,2])
  weight <- adj[which(adj > 0)]
  type  <- match.arg(type)
  .type <- unlist(strsplit(type, "-"))
  x <- switch(.type[1], "out" = out_str, "in" = in_str)[vert_from]
  y <- switch(.type[2], "out" = out_str, "in" = in_str)[vert_to]
  weighted.cor <- function(x, y, w) {
    mean_x <- stats::weighted.mean(x, w)
    mean_y <- stats::weighted.mean(y, w)
    var_x <- sum((x - mean_x)^2 * w)
    var_y <- sum((y - mean_y)^2 * w)
    return(sum(w * (x - mean_x) * (y - mean_y)) / 
             sqrt(var_x * var_y))
  }
  return(weighted.cor(x, y, weight))
}

#' @importFrom wdm wdm
NULL

#' Return four directed, weighted assortativity coefficient with a given
#' network.
#'
#' @param netmat A three column matrix to use for assortativity calculation.
#' @param directed Logical. Whether the edges will be considered as directed. If
#'   FALSE, the input network will be considered as undirected.
#' @param weighted Logical. Whether the edges will be considered as weighted. If
#'   FALSE, values of the third column will be considered as 1.
#'
#' @return Assortativity coefficient for undirected network, or four directed
#'   assortativity coefficients for directed network.
#' @export
#'
#' @examples
#' net <- rpanet(10^3)
#' result <- edge_assort(net, directed = TRUE)
#' net <- rpanet(10^3, directed = FALSE)
#' result <- edge_assort(net, directed = FALSE)
edge_assort <- function(netmat, directed = TRUE, weighted = TRUE) {
  if (! directed) netmat <- rbind(netmat, netmat[, c(2, 1, 3)])
  if (! weighted) netmat[, 3] <- 1
  numnode <- max(netmat[, c(1, 2)])
  outs <- ins <- rep(0, numnode)
  dataf <- data.frame(netmat)
  colnames(dataf) <- c('x', 'y', 'w')
  touts <- aggregate(w ~ x, data = dataf, FUN = 'sum')
  tins <- aggregate(w ~ y, data = dataf, FUN = 'sum')
  outs[touts[, 1]] <- touts[, 2]
  ins[tins[, 1]] <- tins[, 2]
  w <- netmat[, 3]
  require(wdm)
  if (! directed) return(wdm(x = outs[netmat[, 1]], y = outs[netmat[, 2]], weights = w, method = 'pearson'))
  result <- list('out-out' = wdm(x = outs[netmat[, 1]], y = outs[netmat[, 2]],
                                 weights = w, method = 'pearson'), 
                 'out-in' = wdm(x = outs[netmat[, 1]], y = ins[netmat[, 2]], 
                                weights = w, method = 'pearson'), 
                 'in-out' = wdm(x = ins[netmat[, 1]], y = outs[netmat[, 2]],
                                weights = w, method = 'pearson'),
                 'in-in' = wdm(x = ins[netmat[, 1]], y = ins[netmat[, 2]],
                               weights = w, method = 'pearson'))
  return(result)
}

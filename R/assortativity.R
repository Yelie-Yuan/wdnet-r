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
#' "in-in", "out-out" or "in-out"? For undirected networks, consider function
#' \code{uw_assort}
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

#' Degree-based centrality
#'
#' Compute the degree centrality measures of the vertices in a weighted and directed 
#' network represented through its adjacency matrix.
#'
#' @usage
#' degree_c(adj, alpha = 1, mode = "out")
#'
#' @param adj is an adjacency matrix of an weighted and directed network
#' @param alpha is a tuning parameter. The value of alpha must be nonnegative. By convetion, 
#' alpha takes a value from 0 to 1 (default).
#' @param mode which mode to compute: "out" (default) or "in"? For undirected networks, this
#' setting is irrelavent.
#'
#' @return a list of node names and associated degree centrality measures
#'
#' @references
#' \itemize{
#' \item Opsahl, T., Agneessens, F., Skvoretz, J. (2010). Node centrality 
#' in weighted networks: Generalizing degree and shortest paths. 
#' \emph{Social Networks}, \textbf{32}, 245--251.
#' \item Zhang, P., Zhao, J. and Yan, J. (2020+) Centrality measures of 
#' networks with application to world input-output tables
#' }
#'
#' @note 
#' Function \code{degree_c} is an extension of function \code{strength} 
#' in package \code{igraph} and an altertanive of function \code{degree_w} in 
#' package \code{tnet}. Function \code{degree_c} uses adjacency matrix as
#' input.
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x){x*sample(3,1)})
#' adj_ER <- matrix(weight_ER,20,20)
#' system.time(mydegree <- degree_c(adj_ER, alpha = 0.8, mode = "in"))
#' 
#' @export

degree_c <- function(adj, alpha = 1, mode = "out"){
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  else{
    if (isSymmetric(adj) == TRUE){
      warnings("The analyzed network is undirected!")
    }
    deg_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    adj_name <- colnames(adj)
    if (is.null(adj_name) == FALSE){
      deg_c_output <- adj_name
    }
    else{
      deg_c_output[,1] <- c(1:dim(adj)[1])
    }
    colnames(deg_c_output) <- c("name","deg_c")
    adj_deg <- adj
    adj_deg[which(adj_deg > 0)] <- 1
    if (mode == "in"){
      deg_c_output[,2] <- colSums(adj)^alpha + colSums(adj_deg)^(1 - alpha)
    }
    if (mode == "out"){
      deg_c_output[,2] <- rowSums(adj)^alpha + rowSums(adj_deg)^(1 - alpha)
    }
    return(deg_c_output)
  }
}

#' Closeness centrality
#'
#' Compute the closeness centrality measures of the vertices in a weighted and directed 
#' network represented through its adjacency matrix.
#' 
#' @usage
#' closeness_c(dj, alpha = 1, mode = "out", method = "harmonic")
#'
#' @param adj is an adjacency matrix of an weighted and directed network
#' @param alpha is a tuning parameter. The value of alpha must be nonnegative. By convetion, 
#' alpha takes a value from 0 to 1 (default).
#' @param mode which mode to compute: "out" (default) or "in"? For undirected networks, this
#' setting is irrelavent.
#' @param method which method to use: "harmonic" (default) or "standard"?
#' @param distance whether to consider the entries in the adjacency matrix as distances or
#' strong connections. The default setting is \code{FALSE}.
#' 
#' @return a list of node names and associated degree centrality measures
#'
#' @references
#' \itemize{
#' \item Newman, M.E.J. (2003). The structure and function of complex
#' networks. \emph{SIAM review}, \textbf{45}(2), 167--256.
#' \item Opsahl, T., Agneessens, F., Skvoretz, J. (2010). Node centrality 
#' in weighted networks: Generalizing degree and shortest paths. 
#' \emph{Social Networks}, \textbf{32}, 245--251.
#' \item Zhang, P., Zhao, J. and Yan, J. (2020+) Centrality measures of 
#' networks with application to world input-output tables
#' }
#'
#' @note 
#' Function \code{closeness_c} is an extension of function \code{closeness} 
#' in package \code{igraph} and function \code{closeness_w} in 
#' package \code{tnet}. 
#'
#' @examples
#' ## Generate a network according to the Erd\"{o}s-Renyi model of order 20
#' and parameter p = 0.3
#' edge_ER <- rbinom(400,1,0.3)
#' weight_ER <- sapply(edge_ER, function(x){x*sample(3,1)})
#' adj_ER <- matrix(weight_ER,20,20)
#' system.time(myclose <- closeness_c(adj_ER, alpha = 0.8, mode = "out"))
#' 
#' @export

closeness_c <- function(adj, alpha = 1, type = "out", method = "harmonic", distance = FALSE){
  require(igraph)
  if (alpha < 0){
    stop("The tuning parameter alpha must be nonnegative!")
  }
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix must be a square matrix!")
  }
  else{
    if (distance == FALSE){
      adj <- 1/adj
      adj[is.infinite(adj)] <- 0
      adj <- adj^alpha
    } else if (distance == TRUE){
      adj <- adj^alpha
    }
    temp_g <- graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)
    closeness_c_output <- matrix(NA, nrow = dim(adj)[1], ncol = 2)
    closeness_c_output[,1] <- c(1:dim(adj)[1])
    colnames(closeness_c_output) <- c("vertex","closeness")
    if (method == "harmonic"){
      temp_d <- 1/distances(temp_g, mode = type, algorithm = "dijkstra")
      temp_d[temp_d == Inf] <- 0
      if (type == "in"){
        closeness_c_output[,2] <- rowSums(temp_d)
      }
      if (type == "out"){
        closeness_c_output[,2] <- rowSums(temp_d)
      }
    }
    if (method == "standard"){
      temp_d <- distances(temp_g, mode = type, algorithm = "dijkstra")
      temp_d[temp_d == Inf] <- 0
      if (type == "in"){
        closeness_c_output[,2] <- 1/rowSums(temp_d)
      }
      if (type == "out"){
        closeness_c_output[,2] <- 1/rowSums(temp_d)
      }
    }
    return(closeness_c_output)
    options(warn) = -1
  }
}


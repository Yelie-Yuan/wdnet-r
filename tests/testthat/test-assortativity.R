test_that("assortativity coefficient", {
  my_equal_test <- function(adj,
                            directed = TRUE,
                            ret_true) {
    g <- igraph::graph_from_adjacency_matrix(
      adj,
      mode = ifelse(directed, "directed", "undirected"),
      weighted = TRUE
    )
    edgelist <- igraph::as_edgelist(g)
    edgeweight <- igraph::E(g)$weight
    netwk1 <- create_wdnet(
      edgelist = edgelist,
      edgeweight = edgeweight,
      directed = directed
    )
    netwk2 <- create_wdnet(
      adj = adj,
      directed = directed,
      weighted = TRUE
    )

    ret1 <- assortcoef(adj = adj, directed = directed)
    ret2 <- assortcoef(netwk = netwk1)
    ret3 <- assortcoef(netwk = netwk2)
    ret4 <- assortcoef(
      edgelist = edgelist,
      edgeweight = edgeweight,
      directed = directed
    )
    expect_equal(round(unlist(ret1), 2), round(unlist(ret2), 2))
    expect_equal(round(unlist(ret1), 2), round(unlist(ret3), 2))
    expect_equal(round(unlist(ret1), 2), round(unlist(ret4), 2))
    if (!missing(ret_true)) {
      expect_equal(round(unlist(ret1), 2), round(unlist(ret_true), 2))
    }
  }
  # directed, unweighted network
  adj <- matrix(c(1, 0, 1, 1, 1, 0, 1, 0, 0), nrow = 3, ncol = 3)
  my_equal_test(
    adj,
    directed = TRUE,
    ret_true = list(
      "outout" = -0.16666,
      "outin" = -0.408248,
      "inout" = -0.61237,
      "inin" = -0.25
    )
  )
  tmp <- abs(
    assortcoef(adj = adj, directed = TRUE)$outin - 
    igraph::assortativity.degree(
      igraph::graph_from_adjacency_matrix(adj, mode = "directed")
    )
  )
  expect_true(tmp < 1e-5)

  # directed, weighted network
  adj[adj > 0] <- c(
    0.8573071, 1.0371858, 1.0613809,
    1.3421739, 1.0267115
  )
  my_equal_test(
    adj,
    directed = TRUE,
    ret_true = list(
      "outout" = -0.2796024,
      "outin" = -0.3560606,
      "inout" = -0.6660364,
      "inin" = 0.2911594
    )
  )

  # undirected, weighted network
  adj <- pmax(adj, t(adj))
  my_equal_test(
    adj,
    directed = FALSE,
    ret_true = -0.164347
  )

  # undirected, unweighted network
  adj[adj > 0] <- 1
  my_equal_test(
    adj,
    directed = FALSE,
    ret_true = -0.3333333
  )
  tmp <- abs(
    assortcoef(adj = adj, directed = FALSE) - 
    igraph::assortativity.degree(
      igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    )
  )
  expect_true(tmp < 1e-5)

  netwk <- create_wdnet(adj = adj, directed = TRUE)
  tmp <- unlist(
      assortcoef(netwk = netwk)
    ) -
    unlist(
      dw_feature_assort(netwk = netwk, f1 = netwk$node.attr$outs, f2 = netwk$node.attr$ins)
    )
  expect_true(sum(abs(tmp)) < 1e-5)
})

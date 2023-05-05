test_that("centrality coefficient", {
  my_test_equal <- function(ret_degree, ret_close, ret_wpr) {
    ret_degree0 <- data.frame(
      "name" = 1:10,
      "degree" = c(
        3.408225, 4.772597, 5.438694, 5.438694,
        2.741101, 4.869629, 5.341661, 3.556923,
        7.045277, 3.556923
      )
    )
    ret_close0 <- data.frame(
      "name" = 1:10,
      "degree" = c(
        8.803582, 9.393707, 7.781299, 6.982305,
        7.687608, 13.872578, 0.000000, 6.724363,
        6.546484, 13.047310
      )
    )
    ret_wpr0 <- data.frame(
      "name" = 1:10,
      "degree" = c(
        0.05079972, 0.11095956, 0.15453477, 0.09727950,
        0.03403035, 0.19727252, 0.08086971, 0.10003958,
        0.09747361, 0.07674069
      )
    )
    temp <- sum(abs(ret_degree - ret_degree0))
    temp <- temp + sum(abs(ret_close - ret_close0))
    temp <- temp + sum(abs(ret_wpr - ret_wpr0))
    expect_lt(abs(temp), 1e-5)
  }
  set.seed(1234)
  adj <- matrix(rbinom(100, 1, 0.3) * sample(3, 100, replace = TRUE),
    nrow = 10
  )
  degree.control <- list(alpha = 0.8, mode = "in")
  closeness.control <- list(
    alpha = 0.8, mode = "out",
    method = "harmonic", distance = FALSE
  )
  wpr.control <- list(
    gamma = 0.85,
    theta = 0.75,
    prior.info = rep(1 / nrow(adj), nrow(adj))
  )
  my_test_equal(
    ret_degree = centrality(
      adj = adj,
      measure = "degree",
      degree.control = degree.control
    ),
    ret_close = centrality(
      adj = adj,
      measure = "closeness",
      closeness.control = closeness.control
    ),
    ret_wpr = centrality(
      adj = adj,
      measure = "wpr",
      wpr.control = wpr.control
    )
  )

  netwk <- create_wdnet(adj = adj, directed = TRUE, weighted = TRUE)
  my_test_equal(
    ret_degree = centrality(
      netwk = netwk,
      measure = "degree",
      degree.control = degree.control
    ),
    ret_close = centrality(
      netwk = netwk,
      measure = "closeness",
      closeness.control = closeness.control
    ),
    ret_wpr = centrality(
      netwk = netwk,
      measure = "wpr",
      wpr.control = wpr.control
    )
  )
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "directed",
    weighted = TRUE
  )
  edgelist <- igraph::as_edgelist(g)
  edgeweight <- igraph::E(g)$weight
  my_test_equal(
    ret_degree = centrality(
      edgelist = edgelist,
      edgeweight = edgeweight,
      measure = "degree",
      degree.control = degree.control
    ),
    ret_close = centrality(
      edgelist = edgelist,
      edgeweight = edgeweight,
      measure = "closeness",
      closeness.control = closeness.control
    ),
    ret_wpr = centrality(
      edgelist = edgelist,
      edgeweight = edgeweight,
      measure = "wpr",
      wpr.control = wpr.control
    )
  )
  netwk <- create_wdnet(
    edgelist = edgelist,
    edgeweight = edgeweight,
    directed = TRUE
  )
  my_test_equal(
    ret_degree = centrality(
      netwk = netwk,
      measure = "degree",
      degree.control = degree.control
    ),
    ret_close = centrality(
      netwk = netwk,
      measure = "closeness",
      closeness.control = closeness.control
    ),
    ret_wpr = centrality(
      netwk = netwk,
      measure = "wpr",
      wpr.control = wpr.control
    )
  )
})

test_that("clustering coefficient", {
  mycc0 <- list(
    "total" = list(
      "localcc" = c(
        0.0000000, 0.4259259, 0.4000000, 0.3750000,
        0.0000000, 0.1926230, 0.0000000, 0.3888889,
        0.2500000, 0.2604167
      ),
      "numtriangles" = c(
        0.0, 11.5, 8.0, 13.5, 0.0, 23.5,
        0.0, 7.0, 12.0, 12.5
      ),
      "globalcc" = 0.2292854
    ),
    "out" = list(
      "localcc" = c(
        0.0000000, 0.0000000, 0.0000000, 0.5000000,
        0.0000000, 0.1442308, 0.0000000, 1.0000000,
        0.0000000, 0.2666667
      ),
      "globalcc" = 0.1910897,
      "numtriangles" = c(
        0.0, 0.0, 0.0, 1.5, 0.0, 7.5,
        0.0, 3.0, 0.0, 8.0
      )
    ),
    "in" = list(
      "localcc" = c(
        0.0000000, 0.5000000, 0.3333333, 0.5000000, 0.0000000,
        0.1500000, 0.0000000, 0.5000000, 0.1666667, 0.0000000
      ),
      "globalcc" = 0.215,
      "numtriangles" = c(
        0.0, 2.5, 4.0, 6.0, 0.0,
        1.5, 0.0, 1.5, 3.0, 0.0
      )
    ),
    "middle" = list(
      "localcc" = c(
        0.0000000, 0.2222222, 0.5000000, 0.1904762, 0.0000000,
        0.2000000, 0.0000000, 0.0000000, 0.4814815, 0.0000000
      ),
      "globalcc" = 0.159418,
      "numtriangles" = c(
        0.0, 2.0, 2.0, 2.0, 0.0, 6.0,
        0.0, 0.0, 6.5, 0.0
      )
    ),
    "cycle" = list(
      "localcc" = c(
        0.0000000, 0.7777778, 0.5000000, 0.3809524, 0.0000000,
        0.2833333, 0.0000000, 0.4166667, 0.1851852, 0.5000000
      ),
      "globalcc" = 0.3043915,
      "numtriangles" = c(
        0.0, 7.0, 2.0, 4.0, 0.0,
        8.5, 0.0, 2.5, 2.5, 4.5
      )
    )
  )

  my_test_equal <- function(x, y) {
    temp <- sapply(
      seq_along(x), function(i) {
        sapply(
          seq_along(i), function(j) {
            sum(abs(x[[i]][[j]] - y[[i]][[j]]))
          }
        )
      }
    )
    temp <- sum(abs(temp))
    expect_equal(round(temp, 5), 0)
  }
  set.seed(1234)
  adj <- matrix(rbinom(100, 1, 0.3) * sample(3, 100, replace = TRUE),
    nrow = 10
  )
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "directed",
    weighted = TRUE
  )
  edgelist <- igraph::as_edgelist(g)
  edgeweight <- igraph::E(g)$weight

  netwk <- create_wdnet(
    edgelist = edgelist, edgeweight = edgeweight,
    directed = TRUE
  )
  mycc <- clustcoef(netwk = netwk, method = "Clemente")
  my_test_equal(mycc, mycc0)

  netwk <- create_wdnet(
    adj = adj, directed = TRUE, weighted = TRUE
  )
  mycc <- clustcoef(netwk = netwk, method = "Clemente")
  my_test_equal(mycc, mycc0)

  mycc <- clustcoef(adj = adj, method = "Clemente")
  my_test_equal(mycc, mycc0)

  mycc <- clustcoef(
    edgelist = edgelist,
    edgeweight = edgeweight, method = "Clemente"
  )
  my_test_equal(mycc, mycc0)
})

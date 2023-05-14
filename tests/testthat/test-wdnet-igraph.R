test_that("multiplication works", {
  set.seed(123)
  netwk0 <- rpanet(nstep = 1e3)
  g0 <- wdnet_to_igraph(netwk0)
  netwk1 <- igraph_to_wdnet(g0)
  g1 <- wdnet_to_igraph(netwk1)
  netwk2 <- igraph_to_wdnet(g1)
  g2 <- wdnet_to_igraph(netwk2)
  tmp <- c("edgelist", "node.attr", "edge.attr", "weighted", "directed")
  expect_true(all(
    igraph::identical_graphs(g0, g1),
    igraph::identical_graphs(g0, g2),
    identical(netwk0[tmp], netwk1[tmp]),
    identical(netwk0[tmp], netwk2[tmp])))
})

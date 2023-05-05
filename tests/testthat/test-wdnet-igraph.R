test_that("multiplication works", {
  set.seed(123)
  netwk0 <- rpanet(nstep = 1e3)
  g0 <- wdnet_to_igraph(netwk0)
  netwk1 <- igraph_to_wdnet(g0)
  g1 <- wdnet_to_igraph(netwk1)

  expect_true(igraph::identical_graphs(g0, g1))
  expect_true(all(
    identical(netwk0$edgelist, netwk1$edgelist),
    identical(netwk0$node.attr, netwk1$node.attr),
    identical(netwk0$edge.attr, netwk1$edge.attr),
    identical(netwk0$weighted, netwk1$weighted),
    identical(netwk0$directed, netwk1$directed)))
})

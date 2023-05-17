test_that("rewire", {
  set.seed(123)
  my_expect_true <- function(nstep, control, directed) {
    set.seed(123)
    edgelist <- rpanet(
      nstep = nstep,
      control = control,
      initial.network = list(
        edgelist = matrix(1:2, nrow = 1),
        directed = directed)
    )$edgelist
    tmp <- which(edgelist[, 1] > edgelist[, 2])
    edgelist[tmp, ] <- edgelist[tmp, c(2, 1)]
    tmp_order <- order(edgelist[, 1], edgelist[, 2])
    edgelist <- edgelist[tmp_order, ]
    adj <- edgelist_to_adj(edgelist = edgelist, directed = directed)
    netwk1 <- create_wdnet(
      edgelist = edgelist,
      directed = directed)
    netwk2 <- create_wdnet(
      adj = adj,
      weighted = FALSE,
      directed = directed
    )
    stopifnot(
      identical(edgelist, netwk1$edgelist) &&
      identical(netwk1$edgelist, netwk2$edgelist)
    )
    if (directed) {
      target.assortcoef <- list("outout" = -0.2, "outin" = 0.2)
      eta <- get_eta_directed(
        edgelist = edgelist,
        target.assortcoef = target.assortcoef
      )$eta
    } else {
      target.assortcoef <- 0.2
      eta <- get_eta_undirected(
        edgelist = edgelist,
        target.assortcoef = target.assortcoef
      )$eta
    }

    set.seed(123)
    ret1 <- dprewire(
      edgelist = edgelist, directed = directed,
      control = list(iteration = 200), eta = eta
    )
    set.seed(123)
    ret2 <- dprewire(
      adj = adj, directed = directed,
      control = list(iteration = 200), eta = eta
    )
    set.seed(123)
    ret3 <- dprewire(
      netwk1,
      control = list(iteration = 200), eta = eta
    )
    set.seed(123)
    ret4 <- dprewire(
      netwk2,
      control = list(iteration = 200), eta = eta
    )
    ret <- identical(ret1, ret2) &&
      identical(ret1, ret3) &&
      identical(ret1, ret4)
    expect_true(ret)
  }
  nstep <- 1e3
  control <- rpa_control_scenario(
    alpha = 0.3,
    beta = 0.4,
    gamma = 0.3
  )
  my_expect_true(
    nstep = nstep,
    control = control,
    directed = TRUE
  )
  my_expect_true(
    nstep = nstep,
    control = control,
    directed = FALSE
  )
})

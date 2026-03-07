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
    expect_equal(edgelist, netwk1$edgelist)
    expect_equal(netwk1$edgelist, netwk2$edgelist)
    if (directed) {
      target.assortcoef <- list("outout" = -0.2, "outin" = 0.2)
      eta <- get_eta_directed(
        edgelist = edgelist,
        target.assortcoef = target.assortcoef
      )$eta
      
      set.seed(123)
      inin.range1 <- dprewire.range(
        edgelist = edgelist, which.range = 'inin',
        target.assortcoef = target.assortcoef
      )$range
      expect_lt(inin.range1[1], -0.6)
      expect_gt(inin.range1[2], 0.8)
      # print(inin.range1)
      inin.range2 <- dprewire.range(
        edgelist = edgelist, which.range = 'inin'
      )$range
      expect_lt(inin.range2[1], -0.7)
      expect_gt(inin.range2[2], 0.9)
      # print(inin.range2)
      expect_lte(inin.range2[1], inin.range1[1])
      expect_gte(inin.range2[2], inin.range1[2])

      inin.range3 <- dprewire.range(
        edgelist = edgelist, which.range = 'inin',
        target.assortcoef = list("outout" = -0.2, "outin" = c(0.15, 0.25))
      )$range
      # print(inin.range3)
      expect_lt(inin.range3[1], -0.6)
      expect_gt(inin.range3[2], 0.8)
      expect_lte(inin.range3[1], inin.range1[1])
      expect_gte(inin.range3[2], inin.range1[2])
      expect_gte(inin.range3[1], inin.range2[1])
      expect_lte(inin.range3[2], inin.range2[2])
    } else {
      target.assortcoef <- 0.2
      eta <- get_eta_undirected(
        edgelist = edgelist,
        target.assortcoef = target.assortcoef
      )$eta
      
      set.seed(123)
      myrange <- dprewire.range(
        edgelist = edgelist
      )$range
      # print(myrange)
      expect_lt(myrange[1], -0.2)
      expect_gt(myrange[2], 0.5)
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
    expect_equal(ret1, ret2)
    expect_equal(ret1, ret3)
    expect_equal(ret1, ret4)
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

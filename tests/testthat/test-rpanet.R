test_that("rpanet with default preference functions", {
  # sample PA networks
  set.seed(1234)
  nstep <- 1e3
  my_weight_sampler <- function(n) {
    rgamma(n, shape = 5, scale = 0.2)
  }
  for (method in c("linear", "binary", "bag", "bagx")) {
    if (method == "linear" | method == "binary") {
      control <- rpa_control_preference(
        ftype = "default",
        sparams = runif(5, 1, 3),
        tparams = runif(5, 1, 3),
        params = runif(2, 1, 3)
      ) +
        rpa_control_scenario(
          alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1
        ) +
        rpa_control_edgeweight(sampler = my_weight_sampler)
    } else if (method == "bag") {
      control <- rpa_control_preference(
        ftype = "default",
        sparams = c(1, 1, 0, 0, 0.3),
        tparams = c(0, 0, 1, 1, 0.3),
        params = c(1, 0.3)
      ) +
        rpa_control_scenario(
          alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1
        )
    } else {
      control <- rpa_control_preference(
        ftype = "default",
        sparams = c(1, 1, 0, 0, 0.2),
        tparams = c(0, 0, 1, 1, 0.2),
        params = c(1, 0.2)
      ) +
        rpa_control_scenario(
          alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1
        ) +
        rpa_control_edgeweight(sampler = my_weight_sampler)
    }
    initial.network1 <- rpanet(1e3,
      initial.network = list(
        edgelist = matrix(1:2, nrow = 1),
        directed = TRUE
      ),
      control = control
    )
    initial.network2 <- rpanet(1e3,
      initial.network = list(
        edgelist = matrix(1:2, nrow = 1),
        directed = FALSE
      ),
      control = control
    )
    net1 <- rpanet(
      control = control, nstep = nstep, initial.network = initial.network1,
      method = method
    )
    net2 <- rpanet(
      control = control, nstep = nstep, initial.network = initial.network2,
      method = method
    )

    # check node preference
    sparams <- control$preference$sparams
    tparams <- control$preference$tparams
    params <- control$preference$params
    ret1.1 <- range(net1$node.attr$spref -
      (sparams[1] * net1$node.attr$outs^sparams[2] +
        sparams[3] * net1$node.attr$ins^sparams[4] +
        sparams[5]))
    ret1.2 <- range(net1$node.attr$tpref -
      (tparams[1] * net1$node.attr$outs^tparams[2] +
        tparams[3] * net1$node.attr$ins^tparams[4] +
        tparams[5]))
    ret2 <- range(
      net2$node.attr$pref - (net2$node.attr$s^params[1] + params[2])
    )
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "default, diff preference", ret, "\n")
    expect_lt(ret, 1e-5)

    # check node strength
    temp1 <- node_strength_cpp(
      net1$edgelist[, 1],
      net1$edgelist[, 2],
      net1$edge.attr$weight,
      max(net1$edgelist),
      TRUE
    )
    temp2 <- node_strength_cpp(
      net2$edgelist[, 1],
      net2$edgelist[, 2],
      net2$edge.attr$weight,
      max(net2$edgelist),
      TRUE
    )
    ret1.1 <- range(net1$node.attr$outs - temp1$outs)
    ret1.2 <- range(net1$node.attr$ins - temp1$ins)
    ret2 <- range(net2$node.attr$s - (temp2$outs + temp2$ins))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "default, diff strength", ret, "\n")
    expect_lt(ret, 1e-5)
  }
})

test_that("rpanet with customized preference functions", {
  # sample PA networks
  set.seed(12345)
  nstep <- 1e4
  for (method in c("linear", "binary")) {
    control <- rpa_control_preference(
      ftype = "customized",
      spref = "outs + pow(ins, 0.5) + 1",
      tpref = "pow(outs, 0.5) + ins + 1",
      pref = "pow(s, 1.5) + 1"
    ) +
      rpa_control_scenario(
        alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1
      ) +
      rpa_control_edgeweight(
        sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
      )
    initial.network1 <- rpanet(1e3,
      initial.network = list(
        edgelist = matrix(1:2, nrow = 1),
        directed = TRUE
      ),
      control = control
    )
    initial.network2 <- rpanet(1e3,
      initial.network = list(
        edgelist = matrix(1:2, nrow = 1),
        directed = FALSE
      ),
      control = control
    )
    net1 <- rpanet(
      control = control, nstep = nstep, initial.network = initial.network1,
      method = method
    )
    net2 <- rpanet(
      control = control, nstep = nstep, initial.network = initial.network2,
      method = method
    )

    # check node preference
    ret1.1 <- range(net1$node.attr$spref -
      (net1$node.attr$outs + net1$node.attr$ins^0.5 + 1))
    ret1.2 <- range(net1$node.attr$tpref -
      (net1$node.attr$outs^0.5 + net1$node.attr$ins + 1))
    ret2 <- range(net2$node.attr$pref - (net2$node.attr$s^1.5 + 1))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "customized, diff preference", ret, "\n")
    expect_lt(ret, 1e-5)

    # check node strength
    temp1 <- node_strength_cpp(
      net1$edgelist[, 1],
      net1$edgelist[, 2],
      net1$edge.attr$weight,
      max(net1$edgelist),
      TRUE
    )
    temp2 <- node_strength_cpp(
      net2$edgelist[, 1],
      net2$edgelist[, 2],
      net2$edge.attr$weight,
      max(net2$edgelist),
      TRUE
    )
    ret1.1 <- range(net1$node.attr$outs - temp1$outs)
    ret1.2 <- range(net1$node.attr$ins - temp1$ins)
    ret2 <- range(net2$node.attr$s - (temp2$outs + temp2$ins))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "customized, diff strength", ret, "\n")
    expect_lt(ret, 1e-5)
  }
})


test_that("rpanet initial network", {
  set.seed(12345)
  ctr1 <- rpa_control_preference(
    ftype = "customized",
    spref = "outs + pow(ins, 0.5) + 1",
    tpref = "pow(outs, 0.5) + ins + 1",
    pref = "pow(s, 1.5) + 1"
  ) + rpa_control_scenario(
    alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1
  ) + rpa_control_edgeweight(
    sampler = function(n) rgamma(n, shape = 5, scale = 0.2)
  ) + rpa_control_newedge(
    sampler = function(n) rpois(n, lambda = 2) + 1
  ) + rpa_control_reciprocal(
    group.prob = c(0.2, 0.4, 0.4),
    recip.prob = matrix(rep(0.5, 9), nrow = 3)
  )

  netwk1 <- rpanet(1e3,
    initial.network = list(
      edgelist = matrix(1:2, nrow = 1),
      directed = TRUE
    ),
    control = ctr1
  )
  netwk2 <- rpanet(1e3,
    initial.network = netwk1,
    control = netwk1$control
  )
  netwk3 <- rpanet(1e3,
    initial.network = list(
      edgelist = matrix(1:2, nrow = 1),
      directed = FALSE
    ),
    control = ctr1
  )
  netwk4 <- rpanet(1e3,
    initial.network = netwk3,
    control = netwk1$control
  )

  # check initial netwk
  check_initial_network <- function(netwk1, netwk2) {
    nedge <- nrow(netwk1$edgelist)
    nnode <- nrow(netwk1$node.attr)
    netwk1$edge.attr$scenario <- 0
    expect_equal(netwk1$edgelist, netwk2$edgelist[1:nedge, ])
    expect_true(all(netwk2$edge.attr$scenario[1:nedge] == 0))
    expect_equal(netwk1$edge.attr$weight, netwk2$edge.attr$weight[1:nedge])
    expect_equal(netwk1$directed, netwk2$directed)
    expect_equal(netwk1$weighted, netwk2$weighted)
    expect_equal(netwk1$control, netwk2$control)
    expect_equal(netwk1$node.attr$group, netwk2$node.attr$group[1:nnode])
    NULL
  }
  check_initial_network(netwk1, netwk2)
  check_initial_network(netwk3, netwk4)
})


test_that("rpanet scenarios", {
  set.seed(12345)
  ctr <- rpa_control_scenario(
    alpha = 0.1, beta = 0.8, gamma = 0.1, beta.loop = FALSE
  )
  netwk1 <- rpanet(1e4,
    control = ctr,
    initial.network = list(
      edgelist = matrix(1:2, nrow = 1),
      directed = TRUE
    )
  )
  netwk2 <- rpanet(1e4,
    control = ctr,
    initial.network = list(
      edgelist = matrix(1:2, nrow = 1),
      directed = FALSE
    )
  )

  check_scenarios <- function(netwk) {
    alpha <- which(netwk$edge.attr$scenario == 1)
    beta <- which(netwk$edge.attr$scenario == 2)
    gamma <- which(netwk$edge.attr$scenario == 3)
    expect_true(all(netwk$edgelist[alpha, 1] > netwk$edgelist[alpha, 2]))
    expect_true(all(netwk$edgelist[beta, 1] != netwk$edgelist[beta, 2]))
    expect_true(all(netwk$edgelist[gamma, 1] < netwk$edgelist[gamma, 2]))
  }

  check_scenarios(netwk1)
  check_scenarios(netwk2)
})

test_that("rpanet node id", {
  set.seed(12345)
  ctr <- rpa_control_scenario(
    alpha = 0.1, beta = 0.8, gamma = 0.1, beta.loop = FALSE
  )
  initial.network <- list(
    edgelist = matrix(c(1010:1050, 3010:3050), ncol = 2), directed = TRUE
  )
  netwk1 <- rpanet(1e4,
    control = ctr,
    initial.network = initial.network
  )
  expect_equal(
    netwk1$edgelist[seq_len(nrow(initial.network$edgelist)), ],
    initial.network$edgelist
  )
})

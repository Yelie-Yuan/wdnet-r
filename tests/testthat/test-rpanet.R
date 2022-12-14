test_that("Test rpanet with default preference functions", {
  # sample PA networks
  set.seed(1234)
  nstep <- 1e5
  for (method in c("linear", "binary", "bag", "bagx")) {
    if (method == "linear" | method == "binary") {
      control <- rpa_control_preference(ftype = "default",
                                        sparams = runif(5, 1, 3),
                                        tparams = runif(5, 1, 3),
                                        params = runif(2, 1, 3)) +
        rpa_control_scenario(alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1) +
        rpa_control_edgeweight(distribution = rgamma, dparams = list(shape = 5, scale = 0.2))
    } else if (method == "bag") {
      control <- rpa_control_preference(ftype = "default",
                                        sparams = c(1, 1, 0, 0, 0.3),
                                        tparams = c(0, 0, 1, 1, 0.3),
                                        params = c(1, 0.3)) +
        rpa_control_scenario(alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1)
    } else {
      control <- rpa_control_preference(ftype = "default",
                                        sparams = c(1, 1, 0, 0, 0.2),
                                        tparams = c(0, 0, 1, 1, 0.2),
                                        params = c(1, 0.2)) +
        rpa_control_scenario(alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1) +
        rpa_control_edgeweight(distribution = rgamma, dparams = list(shape = 5, scale = 0.2))
    }
    initial.network1 <- rpanet(1e3, directed = TRUE, control = control)
    initial.network2 <- rpanet(1e3, directed = FALSE, control = control)
    net1 <- rpanet(control = control, nstep = nstep, initial.network = initial.network1,
                   directed = TRUE, method = method)
    net2 <- rpanet(control = control, nstep = nstep, initial.network = initial.network2,
                   directed = FALSE, method = method)
    
    # check node preference
    sparams <- control$preference$sparams
    tparams <- control$preference$tparams
    params <- control$preference$params
    ret1.1 <- range(net1$node.attribute$spref - 
                      (sparams[1] * net1$node.attribute$outstrength^sparams[2] + 
                         sparams[3] * net1$node.attribute$instrength^sparams[4] + 
                         sparams[5]))
    ret1.2 <- range(net1$node.attribute$tpref - 
                      (tparams[1] * net1$node.attribute$outstrength^tparams[2] + 
                         tparams[3] * net1$node.attribute$instrength^tparams[4] + 
                         tparams[5]))
    ret2 <- range(net2$node.attribute$pref - (net2$node.attribute$strength^params[1] + params[2]))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "default, diff preference", ret, "\n")
    expect_lt(ret, 1e-5)
    
    # check node strength
    temp1 <- node_strength_cpp(net1$edgelist[, 1], 
                               net1$edgelist[, 2],
                               net1$edgeweight, 
                               max(net1$edgelist), 
                               TRUE)
    temp2 <- node_strength_cpp(net2$edgelist[, 1], 
                               net2$edgelist[, 2],
                               net2$edgeweight, 
                               max(net2$edgelist),
                               TRUE)
    ret1.1 <- range(net1$node.attribute$outstrength - temp1$outstrength)
    ret1.2 <- range(net1$node.attribute$instrength - temp1$instrength)
    ret2 <- range(net2$node.attribute$strength - (temp2$outstrength + temp2$instrength))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "default, diff strength", ret, "\n")
    expect_lt(ret, 1e-5)
  }
})

test_that("Test rpanet with customized preference functions", {
  # sample PA networks
  set.seed(12345)
  nstep <- 1e5
  for (method in c("linear", "binary")) {
    control <- rpa_control_preference(ftype = "customized",
                                      spref = "outs + pow(ins, 0.5) + 1",
                                      tpref = "pow(outs, 0.5) + ins + 1",
                                      pref = "pow(s, 1.5) + 1") +
      rpa_control_scenario(alpha = 0.2, beta = 0.4, gamma = 0.2, xi = 0.1, rho = 0.1) +
      rpa_control_edgeweight(distribution = rgamma, dparams = list(shape = 5, scale = 0.2))
    initial.network1 <- rpanet(1e3, directed = TRUE, control = control)
    initial.network2 <- rpanet(1e3, directed = FALSE, control = control)
    net1 <- rpanet(control = control, nstep = nstep, initial.network = initial.network1,
                   directed = TRUE, method = method)
    net2 <- rpanet(control = control, nstep = nstep, initial.network = initial.network2,
                   directed = FALSE, method = method)
    
    # check node preference
    ret1.1 <- range(net1$node.attribute$spref - 
                      (net1$node.attribute$outstrength + net1$node.attribute$instrength^0.5 + 1))
    ret1.2 <- range(net1$node.attribute$tpref - 
                      (net1$node.attribute$outstrength^0.5 + net1$node.attribute$instrength + 1))
    ret2 <- range(net2$node.attribute$pref - (net2$node.attribute$strength^1.5 + 1))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "customized, diff preference", ret, "\n")
    expect_lt(ret, 1e-5)
    
    # check node strength
    temp1 <- node_strength_cpp(net1$edgelist[, 1],
                               net1$edgelist[, 2],
                               net1$edgeweight,
                               max(net1$edgelist),
                               TRUE)
    temp2 <- node_strength_cpp(net2$edgelist[, 1],
                               net2$edgelist[, 2],
                               net2$edgeweight,
                               max(net2$edgelist),
                               TRUE)
    ret1.1 <- range(net1$node.attribute$outstrength - temp1$outstrength)
    ret1.2 <- range(net1$node.attribute$instrength - temp1$instrength)
    ret2 <- range(net2$node.attribute$strength - (temp2$outstrength + temp2$instrength))
    ret <- max(abs(c(ret1.1, ret1.2, ret2)))
    # cat("\n", "customized, diff strength", ret, "\n")
    expect_lt(ret, 1e-5)
  }
})
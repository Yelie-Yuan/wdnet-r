test_that("rpanet set.seed", {
  load(test_path("testdata", "rpanet-setseed.RData"))

  set.seed(123)

  netwk_binary_new <- rpanet(
    nstep = nstep,
    control = control1,
    method = "binary"
  )

  netwk_linear_new <- rpanet(
    nstep = nstep,
    control = control1,
    method = "linear"
  )

  netwk_bag_new <- rpanet(
    nstep = nstep,
    control = control0,
    method = "bag"
  )

  netwk_bagx_new <- rpanet(
    nstep = nstep,
    control = control0,
    method = "bagx"
  )

  expect_true(
    all(
      identical(netwk_binary_new, netwk_binary),
      identical(netwk_linear_new, netwk_linear),
      identical(netwk_bag_new, netwk_bag),
      identical(netwk_bagx_new, netwk_bagx)
    )
  )
})
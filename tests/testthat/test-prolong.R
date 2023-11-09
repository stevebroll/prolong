test_that("prolong works", {
  data("sim_metabs")
  data("sim_outcome")

  expect_error(prolong(
    sim_metabs,
    sim_outcome,
    lambda2 = 0.25,
    lambdar = 0.15
  ),
  NA)
  expect_error(prolong(
    sim_metabs,
    sim_outcome,
    groups = FALSE,
    lambda2 = 0.25,
    lambdar = 0.15
  ),
  NA)
  expect_error(prolong(
    sim_metabs,
    sim_outcome,
    foldids = rep(rep(1:5, 3), each = 3),
    lambda2 = 0.25,
    lambdar = 0.15
  ),
  NA)
  expect_error(prolong(
    sim_metabs,
    sim_outcome,
    groups = rep(1:ncol(sim_metabs), each = 6),
    lambda2 = 0.25,
    lambdar = 0.15
  ),
  NA)
  expect_error(prolong(
    sim_metabs,
    sim_outcome,
    lambda1 = 0.01,
    lambda2 = 0.25,
    lambdar = 0.15
  ),
  NA)
  expect_error(suppressWarnings(prolong(
    sim_metabs,
    sim_outcome,
    groups = rep(1:ncol(sim_metabs), each = 6),
    optimvals = c(1, 1)
  )),
  NA)
})

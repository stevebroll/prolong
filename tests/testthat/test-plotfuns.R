test_that("delta_heatmap works", {
  data("sim_metabs")
  expect_error(delta_heatmap(sim_metabs), NA)
  expect_error(delta_heatmap(sim_metabs, timediff = "3-2"), NA)
  expect_error(delta_heatmap(sim_metabs, timediff = "4-1"), NA)
  expect_error(delta_heatmap(sim_metabs, timediff = "4-3"), NA)
  expect_error(delta_heatmap(sim_metabs, interactive = FALSE), NA)
  expect_error(delta_heatmap(sim_metabs, grayscale = TRUE), NA)
})

test_that("plot_trajectories works", {
  data("sim_metabs")
  expect_error(plot_trajectories(sim_metabs, selected = 1:20), NA)
  expect_error(plot_trajectories(sim_metabs, selected = 1:20, timelabs = c("0", "2", "4", "6")), NA)
  expect_error(plot_trajectories(sim_metabs, selected = 1:20, colors = "#000000"), NA)
  expect_error(plot_trajectories(sim_metabs, selected = 1:20, colors = grDevices::hcl.colors(nrow(sim_metabs))), NA)

  data("sim_outcome")
  promod <- prolong(sim_metabs, sim_outcome, lambda2 = 0.2208, lambdar = 0.3313)
  expect_error(plot_trajectories(sim_metabs, object = promod), NA)
  expect_error(plot_trajectories(sim_metabs, object = promod, selected = 1:20), NA)
})


test_that("delta_scatter works", {
  data("sim_metabs")
  expect_error(delta_scatter(sim_metabs), NA)
  expect_error(delta_scatter(sim_metabs, timediff2 = "4-1"), NA)
  expect_error(delta_scatter(sim_metabs, timediff1 = "4-1", timediff2 = "3-2"), NA)
  expect_error(delta_scatter(sim_metabs, fisherz = FALSE), NA)
  expect_error(delta_scatter(sim_metabs, interactive = FALSE), NA)
  expect_error(delta_scatter(sim_metabs, digits = 5), NA)

  expect_error(delta_scatter(sim_metabs, timediff1 = "2-1", timediff2 = "3-2", timediff3 = "4-3"), NA)
  expect_error(delta_scatter(sim_metabs, timediff3 = "4-3", fisherz = FALSE), NA)
  expect_error(delta_scatter(sim_metabs, timediff3 = "4-3", interactive = FALSE), NA)
  expect_error(delta_scatter(sim_metabs, timediff3 = "4-3", digits = 5), NA)
})

test_that("delta_network works", {
  data("sim_metabs")
  expect_error(delta_network(sim_metabs), NA)
  expect_error(delta_network(sim_metabs, timediff = "3-2"), NA)
  expect_error(delta_network(sim_metabs, timediff = "3-2", corr_thresh = 0.9), NA)
  expect_error(delta_network(sim_metabs, partial = FALSE), NA)
  expect_error(delta_network(sim_metabs, partial = FALSE, corr_thresh = 0.3), NA)
  expect_error(delta_network(sim_metabs, interactive = FALSE), NA)
  expect_error(delta_network(sim_metabs, partial = FALSE, interactive = FALSE), NA)
  expect_error(delta_network(sim_metabs, partial = FALSE, interactive = FALSE,
                             method = "spearman"), NA)
})

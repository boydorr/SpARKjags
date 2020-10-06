library(SpARKjags)

context("Test plot_density()")

# Load in data
data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

test_that("a runjags object is read in", {
  # Run model and save output to tmp/a.rds

  location <- system.file("test_data", "a.rds", package = "SpARKjags")
  res <- get_model(location)

  expect_silent(g <- plot_density(model = res,
                                  data = data,
                                  var.regex = get_vars(res),
                                  params = get_params(),
                                  labels = get_labels(data)))

  expect_equal(class(g), c("egg", "gtable", "gTree", "grob", "gDesc"))
})

unlink("tmp_test_dir", recursive = TRUE)



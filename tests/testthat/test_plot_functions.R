library(SpARKjags)

context("Test plot_density()")

# Load in data
data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

location <- system.file("test_data", "a.rds", package = "SpARKjags")
res <- get_model(location)

test_that("plot_density runs without error", {
  expect_silent(g <- plot_density(model = res,
                                  data = data))

  expect_equal(class(g), c("egg", "gtable", "gTree", "grob", "gDesc"))
})

test_that("plot_caterpillar runs without error", {
  expect_silent(g <- plot_caterpillar(model = res))

  expect_equal(class(g), c("gg", "ggplot"))
})

test_that("plot_antibiotics runs without error", {
  expect_silent(g <- plot_antibiotics(model = res,
                                      data = data))

  expect_equal(class(g), c("egg", "gtable", "gTree", "grob", "gDesc"))
})



# write cache_model() function to replace run model
# write delete_all_caches() function

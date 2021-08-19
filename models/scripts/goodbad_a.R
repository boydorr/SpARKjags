library(SpARKjags)
library(runjags)

#' Prior values for precision: `tau.class ~ dgamma(0.001, 0.001)`

x <- seq(-0, 10, by = 0.02)
y <- dgamma(x, shape = 0.001, rate = 0.001)
plot(x, y, type = "l", main = "Prior value for precision")

#' Probability of being in the bad group: `prob.of.bad.hosp ~ dbeta(1, 1)`

x <- seq(0, 1, by = 0.02)
y <- dbeta(x, shape1 = 1, shape2 = 1)
plot(x, y, type = "l")

#' Difference between distribution means of the "good" and "bad" effect:
#' `diff ~ dgamma(0.001, 0.001)` and
#' `intercept.plus <- intercept + diff`

x <- seq(-0, 10, by = 0.02)
y <- dgamma(x, shape = 0.001, rate = 0.001)
plot(x, y, type = "l")

#' # Prior value for intercept: `intercept ~ dnorm(0, 0.001)`

x <- seq(-10, 10, by = 0.02)
y <- dnorm(x, mean = 0, sd = 1/0.001)
plot(x, y, type = "l", main = "Prior value for intercept")

#' # Probability of being resistant in the good group (less resistances):
#' `antibiotic.class.effect[a, 1] ~ dnorm(intercept, tau.class)`

x <- seq(-10, 10, by = 0.02)
y <- dnorm(x, mean = 0, sd = 1/0.001)
plot(x, y, type = "l", main = "Prior value for intercept")



data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

directory <- "goodbad_models"

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a.R"),
                            save_to = "test_a",
                            thin = 10)















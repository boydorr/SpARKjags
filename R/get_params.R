#' get_params
#'
#' @export
#'
get_params <- function() {
  list(`probability of being in the bad group` = "prob.of.bad",
       `probability of resistance (good group)` = ",1",
       `probability of resistance (bad group)` = ",2",
       "intercept", "sd")
}

#' get_params2
#'
#' @export
#'
get_params2 <- function() {
  list(`probability of being in the bad group` = "prob.of.bad",
       `probability of resistance` = "$a",
       other = c("intercept", "sd"))
}

good <- function() 1
bad <- function() 2

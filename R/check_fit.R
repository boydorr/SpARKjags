#' check_fit
#'
#' @param model model
#'
#' @export
#'
check_fit <- function(model) {
  # Generate a pointwise log likelihood for each data point, for each
  # mcmc iteration
  tmp <- lapply(as.mcmc.list(model), function(x)
    as.matrix(x[, grep("_like", colnames(x))]))

  # stick chains together
  ploglik <- array(dim = c(dim(tmp[[1]])[1], 2, dim(tmp[[1]])[2]))
  ploglik[,1,] <- tmp[[1]]
  ploglik[,2,] <- tmp[[2]]

  # total number of y's should be 1019
  w <- loo::waic(ploglik)


  # psis loo (a more robust version of waic)
  p <- loo::loo(ploglik, r_eff = loo::relative_eff(exp(ploglik)))

  c(waic = w$estimates['waic', 'Estimate'],
    psis_loo = p$estimates['looic', 'Estimate'])
}

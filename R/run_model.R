#' run_model
#'
#' @export
#'
run_model <- function(data, model_name, directory, thin) {
  if(missing(thin)) {
    results <- run.jags(paste0(directory, "/", model_name, ".R"),
                        data = data$jags,
                        n.chains = 2,
                        silent.jags = T)
  } else {
    results <- run.jags(paste0(directory, "/", model_name, ".R"),
                        data = data$jags,
                        n.chains = 2,
                        silent.jags = T,
                        thin = thin)
  }

  if(!results$summary.available) {
    results <- add.summary(results)
  }

  saveRDS(results, paste0(directory, "/", model_name, ".rds"))
}
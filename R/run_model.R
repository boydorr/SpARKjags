#' run_model
#'
#' @param data data
#' @param model_name model_name
#' @param directory directory
#' @param thin thin
#'
#' @export
#'
run_model <- function(data, model_name, directory, thin) {
  # Find model
  location <- system.file(directory, paste0(model_name, ".R"),
                          package = "SpARKcarbapenem")

  # Run model
  if(missing(thin)) {
    results <- runjags::run.jags(location,
                        data = data$jags,
                        n.chains = 2,
                        silent.jags = T)
  } else {
    results <- runjags::run.jags(location,
                        data = data$jags,
                        n.chains = 2,
                        silent.jags = T,
                        thin = thin)
  }

  if(!results$summary.available) {
    results <- runjags::add.summary(results)
  }

  # Save results
  save_here <- file.path(dirname(location), paste0(model_name, ".rds"))
  saveRDS(results, save_here)
}

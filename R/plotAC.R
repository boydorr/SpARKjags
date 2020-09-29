#' plotAC
#'
#' @export
#'
plotAC <- function(model,
                   var.regex,
                   filename) {

  mcmc.model <- coda::as.mcmc.list(model)
  ind <- grepl(var.regex, colnames(mcmc.model[[1]]))
  mcmc.model <- mcmc.model[, ind, drop = F]

  autocorr.res <- lapply(seq_along(mcmc.model), function(x) {
    tmp <- coda::autocorr.diag(mcmc.model[[x]], lags = 0:50) %>%
      data.frame() %>%
      tibble::rownames_to_column("lag") %>%
      mutate(chain = x)
  }) %>%
    do.call(rbind.data.frame, .) %>%
    reshape2::melt(id.vars = c("lag", "chain"),
                   variable.name = "parameter") %>%
    dplyr::mutate(lag = gsub("Lag ", "", lag),
                  lag = as.numeric(lag))

  g <- autocorr.res %>% ggplot2::ggplot() + ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~parameter) +
    ggplot2::geom_line(ggplot2::aes(x = lag, y = value,
                                    group = chain, colour = chain)) +
    ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed",
                        colour = "red", alpha = 0.6) +
    ggplot2::labs(x = "Lag", y = "Autocorrelation") +
    ggplot2::theme(legend.position = "none")

  if(missing(filename)) {
    g
  } else {
    n <- length(unique(tmp$Parameter))
    ggplot2::ggsave(filename, ., height = 3*n, width = 6, limitsize = F)
  }
}


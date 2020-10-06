#' plot_caterpillar
#'
#' @param model model
#' @param var.regex var.regex
#' @param filename filename
#'
#' @export
#'
plot_caterpillar <- function(model,
                             var.regex,
                             filename) {

  model.ggs <- model %>%
    coda::as.mcmc.list() %>%
    ggmcmc::ggs(family = var.regex) %>%
    dplyr::mutate(Chain = as.factor(.data$Chain))

  g <- ggplot2::ggplot(model.ggs) + ggplot2::theme_minimal() +
    ggplot2::geom_line(ggplot2::aes_string(x = "Iteration", y = "value",
                                           group = "Chain", colour = "Chain"),
                       alpha = 0.8) +
    ggplot2::facet_wrap(~ Parameter, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("grey", "#2c7fb8")) +
    ggplot2::labs(title = "Trace plot(s)") +
    ggplot2::theme(legend.position = "none")

  if(missing(filename)) {
    g
  } else {
    n <- length(unique(model.ggs$Parameter))
    ggplot2::ggsave(filename, g, height = n/ncol, limitsize = F)
  }
}

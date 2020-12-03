#' Potential Scale Reduction Factor
#'
#' The Gelman and Rubin (1992) Potential Scale Reduction Factor (PSRF, or
#' \eqn{\hat{R}}) is used to assess MCMC convergence. If the chains have not
#' converged to a common distribution, the \eqn{\hat{R}} statistic will be
#' greater than one. Values close to 1 are considered good.
#'
#' @param model a \code{runjags} object containing model results
#' @param plot a \code{boolean} object specifying whether or not a plot should
#' be generated
#'
#' @export
#'
testPSRF <- function(model, plot = T) {
  tmp <- model$psrf$psrf %>% data.frame()
  tmp <- cbind.data.frame(Parameter = rownames(tmp), tmp) %>%
    dplyr::filter(!grepl("^bad.", .data$Parameter )) %>%
    mutate(test = case_when(.data$Point.est. >= 1.1 ~ "bad",
                            T ~ "good"))
  if(plot) {
    # plotly::plot_ly(tmp, x = ~Point.est., y = ~Parameter) %>%
    #   plotly::add_trace(type = "scatter", mode = "markers") %>%
    #   plotly::layout(title = "Potential Scale Reduction Factors",
    #                  xaxis = list(title = "PSRF"),
    #                  yaxis = list(title = "Parameter"),
    #                  shapes = list(type = "line",
    #                                y0 = 0,
    #                                y1 = 1,
    #                                yref = "paper",
    #                                x0 = 1.1,
    #                                x1 = 1.1,
    #                                line = list(color = "grey",
    #                                            dash = "dot"))) %>%
    #   return()

    ggplot(tmp) + theme_minimal() +
      geom_point(aes_string(x = "Point.est.", y = "Parameter",
                            colour = "test")) +
      geom_vline(xintercept = 1.1, linetype = "dashed") +
      labs(x = bquote(hat(R)), title = "Potential Scale Reduction Factors") +
      theme(legend.position = "none") %>%
      return()
  } else {
    tmp <- tmp %>%
      filter(.data$test == "bad") %>%
      select(-.data$test)
    if(nrow(tmp) == 0)
      message("All PSRF values < 1.1") else
        return(tmp)
  }
}

#' Effective Sample Size
#'
#' The effective sample size (ESS) is an estimate of the sample size required
#' to achieve the same level of precision if that sample was a simple random
#' sample. In Bayesian statistics, it is common to use the posterior draws from
#' Markov chain Monte Carlo (MCMC) for statistical inference. The MCMC process
#' causes the draws to be correlated. This means that the effective sample size
#' is generally lower than the number of draws. For this reason, the effective
#' sample size – rather than the actual sample size – is typically used when
#' determining if an MCMC model has converged.
#'
#' @param model a \code{runjags} object containing model results
#' @param plot a \code{boolean} object specifying whether or not a plot should
#' be generated
#'
#' @return If any SSEF values are less than or equal to 300, these variables
#' will be returned in a \code{data.frame} along with their values. If all
#' values are greater than 300, a message will state "All SSeff values > 300".
#' @export
#'
testSSEF <- function(model, plot = T) {
  # Label parameters with SSeff less than or equal to 300 as "bad"
  tmp <- data.frame(SSeff = model$mcse$sseff)
  tmp <- cbind.data.frame(tmp, Parameter = rownames(tmp)) %>%
    dplyr::filter(!grepl("^bad.", Parameter)) %>%
    dplyr::mutate(test = dplyr::case_when(.data$SSeff <= 300 ~ "bad",
                                          T ~ "good"))

  if(plot) {
    # plotly::plot_ly(tmp, x = ~SSeff, y = ~Parameter) %>%
    #   plotly::add_trace(type = "scatter", mode = "markers") %>%
    #   plotly::layout(title = "Effective sample size",
    #                  xaxis = list(title = "SSeff"),
    #                  yaxis = list(title = "Parameter"),
    #                  shapes = list(type = "line",
    #                                y0 = 0,
    #                                y1 = 1,
    #                                yref = "paper",
    #                                x0 = 300,
    #                                x1 = 300,
    #                                line = list(color = "grey",
    #                                            dash = "dot"))) %>%
    #   return()

    ggplot2::ggplot(tmp) + ggplot2::theme_minimal() +
      ggplot2::geom_point(ggplot2::aes_string(x = "SSeff", y = "Parameter",
                                              colour = "test")) +
      ggplot2::geom_vline(xintercept = 300, linetype = "dashed") +
      ggplot2::labs(x = "SSeff", title = "Effective sample size") +
      ggplot2::theme(legend.position = "none") %>%
      return()
  } else {
    tmp <- tmp %>%
      dplyr::filter(.data$test == "bad") %>%
      dplyr::select(-.data$test)
    if(nrow(tmp) == 0)
      message("All SSeff values > 300") else
        return(tmp)
  }
}


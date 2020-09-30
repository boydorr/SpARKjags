#' testSSEF
#'
#' @param model model
#' @param plot plot
#'
#' @export
#'
testSSEF <- function(model, plot = F) {
  tmp <- model$mcse$sseff %>%
    data.frame(SSeff = .data)

  tmp <- cbind.data.frame(tmp, Parameter = rownames(tmp)) %>%
    mutate(test = case_when(.data$SSeff <= 300 ~ "bad",
                            T ~ "good"))
  if(plot) {
    ggplot(tmp) + theme_minimal() +
      geom_point(aes_string(x = "SSeff", y = "Parameter", colour = "test")) +
      geom_vline(xintercept = 300, linetype = "dashed") +
      labs(x = "SSeff", title = "Effective sample size") +
      theme(legend.position = "none") %>%
      return()
  } else {
    tmp <- tmp %>%
      filter(.data$test == "bad") %>%
      select(-.data$test)
    if(nrow(tmp) == 0)
      message("All SSeff values > 300") else
        return(tmp)
  }
}

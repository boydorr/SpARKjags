#' testSSEF
#'
#' @export
#'
testSSEF <- function(model, plot = F) {
  tmp <- model$mcse$sseff %>%
    data.frame(SSeff = .) %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(test = case_when(SSeff <= 300 ~ "bad",
                            T ~ "good"))
  if(plot) {
    ggplot(tmp) + theme_minimal() +
      geom_point(aes(x = SSeff, y = Parameter, colour = test)) +
      geom_vline(xintercept = 300, linetype = "dashed") +
      labs(x = "SSeff", title = "Effective sample size") +
      theme(legend.position = "none") %>%
      return()
  } else {
    tmp %<>%
      filter(test == "bad") %>%
      select(-test)
    if(nrow(tmp) == 0)
      message("All SSeff values > 300") else
        return(tmp)
  }
}

#' testPSRF
#'
#' @param model model
#' @param plot plot
#'
#' @export
#'
testPSRF <- function(model, plot = F) {
  tmp <- model$psrf$psrf %>% data.frame()
  tmp <- cbind.data.frame(Parameter = rownames(tmp), tmp) %>%
    dplyr::filter(!grepl("^bad.", .data$Parameter )) %>%
    mutate(test = case_when(.data$Point.est. >= 1.1 ~ "bad",
                            T ~ "good"))
  if(plot)
  {
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

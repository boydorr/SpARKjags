#' testPSRF
#'
#' @export
#'
testPSRF <- function(model, plot = F) {
  tmp <- model$psrf$psrf %>%
    data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    dplyr::filter(!grepl("^bad.", Parameter )) %>%
    mutate(test = case_when(Point.est. >= 1.1 ~ "bad",
                            T ~ "good"))
  if(plot)
  {
    ggplot(tmp) + theme_minimal() +
      geom_point(aes(x = Point.est., y = Parameter, colour = test)) +
      geom_vline(xintercept = 1.1, linetype = "dashed") +
      labs(x = bquote(hat(R)), title = "Potential Scale Reduction Factors") +
      theme(legend.position = "none") %>%
      return()
  } else {
    tmp %<>%
      filter(test == "bad") %>%
      select(-test)
    if(nrow(tmp) == 0)
      message("All PSRF values < 1.1") else
        return(tmp)
  }
}
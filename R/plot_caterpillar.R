#' Generate caterpillar plot
#'
#' Generates a caterpillar plot (trace plot). If a file already exists, it is
#' read into memory.
#'
#' @param model a \code{runjags} object containing model results
#' @param save_to a \code{string} specifying the location in which the output
#' (an \code{rds} file containing the density plot) should be saved
#' @param filename a \code{string} specifying the name of the output (an
#' \code{rds} file containing the density plot)
#' @param save_pdf a \code{boolean} specifying whether or not a pdf should be
#' saved
#'
#' @return Returns a caterpillar plot and saves an rds file
#' @export
#'
plot_caterpillar <- function(model,
                             save_to,
                             filename,
                             save_pdf = FALSE) {
  # Is the plot cached?
  filepath <- file.path(save_to, filename)

  if(file.exists(filepath)) {
    output <- readRDS(filepath)

  } else {
    var.regex <- get_vars(model)
    model.ggs <- model %>%
      coda::as.mcmc.list() %>%
      ggmcmc::ggs(family = var.regex) %>%
      dplyr::mutate(Chain = as.factor(.data$Chain))

    output <- ggplot2::ggplot(model.ggs) + ggplot2::theme_minimal() +
      ggplot2::geom_line(ggplot2::aes_string(x = "Iteration", y = "value",
                                             group = "Chain", colour = "Chain"),
                         alpha = 0.8) +
      ggplot2::facet_wrap(~ Parameter, scales = "free_y") +
      ggplot2::scale_color_manual(values = c("grey", "#2c7fb8")) +
      ggplot2::labs(title = "Trace plot(s)") +
      ggplot2::theme(legend.position = "none")

    # Save (cache) rds file
    if(!file.exists(save_to)) dir.create(save_to, recursive = TRUE)
    saveRDS(output, filepath)
  }

  # Save pdf
  if(save_pdf) {
    n <- length(unique(model.ggs$Parameter))
    pdf_filename <- gsub("rds", "pdf", filepath)
    ggplot2::ggsave(pdf_filename, output, height = n/ncol, limitsize = F)
  }

  output
}

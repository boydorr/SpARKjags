#' plot_jags
#'
#' @export
#'
plot_jags <- function(res, var, labels, intercept, condense = F) {

  if(condense) {

    if(is.null(res$summaries))
      res %<>% add.summary()

    d1 <- res %>%
      coda::as.mcmc.list() %>%
      ggmcmc::ggs() %>%
      dplyr::filter(Parameter %in% var) %>%
      dplyr::group_by(Parameter) %>%
      dplyr::summarise(ymin = min(value),
                       lower = quantile(value, probs = .25),
                       middle = mean(value),
                       upper = quantile(value, probs = .75),
                       ymax = max(value))

    if(!missing(labels))
      d1 %<>% merge(labels) %>%
      dplyr::mutate(Label = factor(Label, levels = labels$Label))

    g1 <- d1 %>%
      ggplot2::ggplot() + ggplot2::coord_flip() + ggplot2::theme_minimal() +
      ggplot2::geom_boxplot(ggplot2::aes(x = Label, ymin = ymin, lower = lower,
                                         middle = middle, upper = upper,
                                         ymax = ymax), stat = "identity") +
      ggplot2::labs(x = "log-odds")

    if(!missing(intercept))
      g1 <- g1 + ggplot2::geom_hline(yintercept = intercept, linetype = "dashed")

    d2 <- res %>%
      coda::as.mcmc.list() %>%
      ggmcmc::ggs() %>%
      dplyr::filter(Parameter %in% var)

    if(!missing(labels))
      d2 %<>% merge(labels) %>%
      dplyr::mutate(Label = factor(Label, levels = labels$Label))

    g2 <- d2 %>%
      ggmcmc::ggs_caterpillar() + ggplot2::theme_minimal()


    if(!missing(intercept))
      g2 <- g2 + ggplot2::geom_vline(xintercept = intercept, linetype = "dashed")

  } else {

    df <- res %>%
      coda::as.mcmc.list() %>%
      ggmcmc::ggs() %>%
      dplyr::filter(Parameter %in% var)

    if(missing(labels)) {
      df %<>% mutate(Label = Parameter)
    } else {
      df %<>% merge(labels) %>%
        dplyr::mutate(Label = factor(Label, levels = labels$Label))
    }

    # Density plot
    g1 <- df %>%
      ggplot2::ggplot() + theme_minimal() +
      ggplot2::geom_density(aes(x = value, group = Chain, fill = Label),
                            alpha = .4) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "log-odds")

    if(missing(labels))
      g1 <- g1 + ggplot2::facet_wrap(~Label, ncol = 1, scales = "free") else
        g1 <- g1 + ggplot2::facet_wrap(~Label, ncol = 1)

    if(!missing(intercept))
      g1 <- g1 + ggplot2::geom_vline(xintercept = intercept, linetype = "dashed")

    # Caterpillar
    g2 <- df %>%
      dplyr::mutate(Chain = as.factor(Chain)) %>%
      ggplot2::ggplot() + theme_minimal() +
      ggplot2::facet_grid(Label ~ ., scale  = 'free_y', switch = 'y') +
      ggplot2::geom_line(ggplot2::aes(x = Iteration, y = value, col = Chain)) +
      ggplot2::theme(legend.position = "none")

    if(!missing(intercept))
      g2 <- g2 + ggplot2::geom_hline(yintercept = intercept, linetype = "dashed")

  }

  cowplot::plot_grid(g1, g2, nrow = 1, rel_widths = c(.6, .4))
}

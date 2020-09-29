#' densityplot
#'
#' @export
#'
densityplot <- function(model,
                        data,
                        var.regex,
                        params,
                        labels = rep(NA, length(params))) {

  df <- import_data(model, data)

  model.ggs <- model %>%
    coda::as.mcmc.list() %>%
    ggmcmc::ggs(family = var.regex)

  # Parameters in model
  parameters <- model.ggs$Parameter %>%
    unique() %>%
    as.character()

  # Relative heights of facets
  rel_heights <- lapply(params, function(this_set) {
    sapply(this_set, function(x)
      parameters[grepl(x, parameters)]) %>%
      unlist() %>%
      unique() %>%
      length()
  }) %>%
    unlist()


  # Plot each facet ---------------------------------------------------------

  plots <- lapply(seq_along(params), function(this_set) {
    # Subset parameters
    plot_these <- sapply(params[this_set], function(x)
      parameters[grepl(x, parameters)]) %>%
      unlist() %>%
      unique()

    title <- names(params)[[this_set]]
    good <- if_else(grepl("(good group)", title), T, F)
    bad <- if_else(grepl("\\(bad group)", title), T, F)

    plot_this <- model.ggs %>%
      dplyr::filter(Parameter %in% plot_these)

    g <- plot_this %>%
      ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::geom_violin(ggplot2::aes(x = Parameter,
                                        y = value),
                           colour = NA,
                           fill = "grey66",
                           scale = "width") +
      ggplot2::stat_summary(ggplot2::aes(x = Parameter, y = value),
                            fun = median,
                            geom = "point",
                            position = ggplot2::position_dodge(width = 1),
                            colour = "grey33") +
      ggplot2::stat_summary(ggplot2::aes(x = Parameter, y = value),
                            fun = function(z) { quantile(z, c(0.25, 0.75)) },
                            geom = "line",
                            colour = "grey33") +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title = title)



    if(title == "probability of being in the bad group") {
      labs <- labels[[1]] %>%
        dplyr::filter(index %in% unique(plot_this$Parameter))
      g <- g + ggplot2::scale_x_discrete(breaks = labs$index,
                                         labels = labs$category)
    }


    # Plot prob of being in the badgroup / prob of resistance --------------

    if(grepl("^probability of resistance", title)) {
      all_labs <- unique(df$label)

      # Extract data
      if(bad) df %<>% dplyr::filter(badgroup == 1)
      if(good) df %<>% dplyr::filter(badgroup == 0)

      # Observed resistance probability
      n <- df %>%
        dplyr::group_by(label) %>%
        dplyr::summarise(n = n())

      if(good | bad) { # goodbad models
        columns <- df %>%
          select(-GUID, -name, -hospital, -clinical,
                 -mean.p.bad, -badgroup, -label) %>%
          colnames()
      } else { # res.a_naive
        columns <- df %>%
          select(-GUID, -name, -hospital, -clinical, -label) %>%
          colnames()
      }

      s <- df %>%
        dplyr::group_by(label) %>%
        dplyr::summarise_at(vars(contains(columns)), sum, na.rm = T)

      labs <- c("Hospital\n(Carriage)",
                "Hospital\n(Clinical)",
                "GP",
                "Outpatients",
                "Volunteers")
      labs <- c(labs, all_labs[!all_labs %in% labs])

      this_many <- length(labs)

      if(this_many <= 12) {
        manual_colours <- setNames(ggthemes::ptol_pal()(this_many),
                                   labs)

      } else {
        colours <- ggthemes::ptol_pal()(12)
        manual_colours <- setNames(colorRampPalette(colours)(this_many),
                                   labs)
      }

      manual_shapes <- setNames(c(16, 16, 1:(this_many-2)), labs)

      labs <- labs[which(labs %in% s$label)]

      probabilities <- merge(n, s) %>%
        dplyr::mutate_at(vars(contains(columns)), ~ . / n) %>%
        reshape2::melt(id.var = "label",
                       measure.vars = columns,
                       variable.name = "classification",
                       value.name = "Probability") %>%
        merge(labels[[this_set]]) %>%
        dplyr::mutate(name = dplyr::case_when(
          grepl("Hospital", label) ~ "Hospital",
          T ~ label),
          label = factor(label),
          label = forcats::fct_relevel(label, labs))

      labs <- levels(probabilities$label)

      manual_colours <- manual_colours[labs]
      manual_shapes <- manual_shapes[labs]

      # Extract labels
      xaxis <- labels[[this_set]] %>%
        tibble::column_to_rownames("index")

      g <- g + ggplot2::geom_point(
        ggplot2::aes(x = index,
                     y = Probability,
                     colour = label,
                     shape = label),
        probabilities) +
        ggplot2::scale_colour_manual(name = "",
                                     values = manual_colours,
                                     labels = labs) +
        ggplot2::scale_shape_manual(name = "",
                                    values = manual_shapes,
                                    labels = labs) +
        ggplot2::scale_x_discrete(breaks = plot_these,
                                  labels = xaxis[plot_these,]) +
        ggplot2::theme(legend.position = "bottom")
      # +
      #   ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))

    }

    # Set axis limits
    if(all(grepl("prob", plot_these))) {
      g <- g + ggplot2::coord_flip(ylim = c(0, 1))
    } else {
      g <- g + ggplot2::coord_flip()
    }
  })

  egg::ggarrange(plots = plots, heights = rel_heights)

}


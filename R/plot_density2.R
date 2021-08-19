#' Generate density plot
#'
#' Generates a violin plot where violins represent the posterior probability of
#' resistance and points represent the calculated probability of resistance to
#' each antibiotic class. If a file already exists, it is read into memory.
#'
#' @param model a \code{runjags} object containing model results
#' @param data a \code{list} containing the data input for the runjags model,
#' generated with \code{jags_data()}
#' @param save_to a \code{string} specifying the location in which the output
#' (an \code{rds} file containing the density plot) should be saved
#' @param filename a \code{string} specifying the name of the output (an
#' \code{rds} file containing the density plot)
#' @param save_pdf a \code{boolean} specifying whether or not a pdf should be
#' saved
#'
#' @return Returns a density plot and saves an rds file
#' @export
#'
plot_density2 <- function(model,
                          data,
                          save_to,
                          filename,
                          save_pdf = FALSE) {

  filepath <- file.path(save_to, filename)

  # Is the plot cached?
  if (file.exists(filepath)) {
    output <- readRDS(filepath)

  } else {
    # data.frame listing SpARK samples, and their resistances to each antibiotic
    # class, as well as the posterior probability of being in the bad group
    # (mean.p.bad), which defines the badgroup (1 if mean.p.bad > 0.5)

    df <- import_data(model, data)

    # Labels (posterior plots are split into good and bad groups as well as
    # clinical and carriage)
    labels <- get_labels(data)
    params <- get_params2()
    var.regex <- get_vars(model)

    model.ggs <- model %>%
      coda::as.mcmc.list() %>%
      ggmcmc::ggs(family = var.regex)

    # Parameters in model
    parameters <- model.ggs$Parameter %>%
      unique() %>%
      as.character()

    # Relative heights of facets
    rel_heights <- c(length(parameters[grepl("^prob", parameters)]),
                     nrow(labels$good),
                     length(parameters[!grepl("prob", parameters)]))

    p1 <- plot1(model.ggs, parameters, labels)
    p2 <- plot2(model, data, parameters, labels, var.regex)
    p3 <- plot3(model.ggs, parameters)
    plots <- list(p1, p2, p3)

    output <- egg::ggarrange(plots = plots, heights = rel_heights)

    # Save (cache) rds file
    if(!file.exists(save_to)) dir.create(save_to, recursive = TRUE)
    saveRDS(output, filepath)
  }

  # Save pdf
  if(save_pdf) {
    pdf_filename <- gsub("rds", "pdf", filepath)
    ggplot2::ggsave(pdf_filename, output, width = 10, height = 10)
  }

  output
}

# -------------------------------------------------------------------------


plot1 <- function(model.ggs, parameters, labels) {
  plot_these <- parameters[grepl("prob.of.bad", parameters)] %>%
    unlist() %>%
    unique() %>%
    c()

  # Get posterior probabilities of resistance to each antibiotic class
  plot_this <- model.ggs %>%
    dplyr::filter(.data$Parameter %in% plot_these)

  yaxis <- labels[[1]] %>%
    dplyr::filter(.data$index %in% unique(plot_this$Parameter))

  # Generate plot
  g <- plot_this %>%
    ggplot2::ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::geom_violin(ggplot2::aes(x = .data$Parameter,
                                      y = .data$value),
                         trim = TRUE,
                         colour = "grey70",
                         fill = "grey90",
                         scale = "width") +
    ggplot2::stat_summary(ggplot2::aes(x = .data$Parameter, y = .data$value),
                          fun = stats::median,
                          geom = "point",
                          position = ggplot2::position_dodge(width = 1),
                          colour = "grey33") +
    ggplot2::stat_summary(ggplot2::aes(x = .data$Parameter, y = .data$value),
                          fun = function(z) {
                            stats::quantile(z, c(0.25, 0.75))
                          },
                          geom = "line",
                          colour = "grey33") +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   # strip.text = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title = "Plot 1",
                  y = "Probability of being in the bad group") +
    ggplot2::scale_x_discrete(breaks = yaxis$index,
                              labels = yaxis$category) +
    ggplot2::coord_flip(ylim = c(0, 1))
}


# -------------------------------------------------------------------------



plot2 <- function(model, data, parameters, labels, var.regex) {

  model.ggs <- model %>%
    coda::as.mcmc.list() %>%
    ggmcmc::ggs(family = var.regex)

  plot_these <- parameters[grepl("^a", parameters)] %>%
    unlist() %>%
    unique() %>%
    c()

  plot_this <- model.ggs %>%
    dplyr::filter(.data$Parameter %in% plot_these) %>%
    dplyr::mutate(Parameter = as.character(.data$Parameter))

  index_clinical <- clinical(data, "index")
  index_carriage <- carriage(data, "index")

  index_good <- good()
  index_bad <- bad()

  if (any(grepl("\\.gp\\.", plot_this$Parameter))) {
    gp <- if_else(clinical(data, "gp"), index_clinical, index_carriage)
    out <- if_else(clinical(data, "out"), index_clinical, index_carriage)
    vol <- if_else(clinical(data, "vol"), index_clinical, index_carriage)

    dat <- plot_this %>%
      # Add 3rd index to denote clinical state for gp, volunteer, and outpatient
      # samples
      dplyr::mutate(Parameter = dplyr::case_when(
        grepl("\\.gp\\.", .data$Parameter) ~ gsub("\\]", paste0(",", gp, "]"),
                                                  .data$Parameter),
        grepl("\\.v\\.", .data$Parameter) ~ gsub("\\]", paste0(",", vol, "]"),
                                                 .data$Parameter),
        grepl("\\.o\\.", .data$Parameter) ~ gsub("\\]", paste0(",", out, "]"),
                                                 .data$Parameter),
        T ~ .data$Parameter)) %>%
      # Extract indices out into new columns corresponding to antimicrobial
      # class (index), goodbad, and clinical state
      dplyr::mutate(index = gsub(".*\\[([0-9]+),[0-9],[0-9]\\]", "\\1",
                                 .data$Parameter),
                    goodbad = gsub(".*\\[[0-9]+,([0-9]),[0-9]\\]", "\\1",
                                   .data$Parameter),
                    clinical = gsub(".*\\[[0-9]+,[0-9],([0-9])\\]", "\\1",
                                    .data$Parameter)) %>%
      dplyr::mutate(index = as.numeric(.data$index),
                    goodbad = if_else(.data$goodbad == index_good, "Good group",
                                      "Bad group"),
                    clinical = if_else(.data$clinical == index_clinical,
                                       "Clinical", "Carriage")) %>%
      dplyr::mutate(Parameter = dplyr::case_when(
        grepl("\\.gp\\.", .data$Parameter) ~ "gp",
        grepl("\\.v\\.", .data$Parameter) ~ "vol",
        grepl("\\.o\\.", .data$Parameter) ~ "out",
        grepl("^ac\\.", .data$Parameter) ~ "hosp"))

    by_clinical <- TRUE

  } else {
    # Extract indices out into new columns corresponding to antimicrobial
    # class (index) and goodbad
    dat <- plot_this %>%
      dplyr::mutate(index = gsub(".*\\[([0-9]+),[0-9]\\]", "\\1",
                                 .data$Parameter),
                    goodbad = gsub(".*\\[[0-9]+,([0-9])\\]", "\\1",
                                   .data$Parameter)) %>%
      dplyr::mutate(index = as.numeric(.data$index),
                    goodbad = dplyr::case_when(
                      .data$goodbad == index_good ~ "Good group",
                      .data$goodbad == index_bad ~ "Bad group"))
  }

  dat <- dat %>%
    # Merge human readable antibiotic class names by index
    left_join(labels$good, by = "index") %>%
    dplyr::select(-.data$index) %>%
    # values are the same across groupings (for ac.prob and gr.prob, for
    # example), so find unique selections
    dplyr::select(-.data$Parameter) %>%
    unique() %>%
    dplyr::mutate(goodbad = factor(.data$goodbad,
                                   levels = c("Good group", "Bad group")))

  # Initialise variables for plotting
  df <- import_data(model, data)

  all_labs <- unique(df$label)

  labs <- c("Hospital\n(Carriage)",
            "Hospital\n(Clinical)",
            "GP",
            "Outpatients",
            "Volunteers")
  labs <- c(labs, all_labs[!all_labs %in% labs])

  this_many <- length(labs)


  # Calculate the observed resistance probability

  # Define which columns should be used when the input is a goodbad model
  # or naive (e.g. res.a_naive)

  columns <- df %>%
    dplyr::select(-.data$GUID, -.data$name, -.data$hospital, -.data$clinical,
                  -.data$mean.p.bad, -.data$badgroup, -.data$label) %>%
    colnames()


  # Count of the number of R and S samples in each hospital-clinical group
  # for each antibiotic class (ignore NAs)

  if (by_clinical) {
    n <- df %>%
      dplyr::group_by(.data$label, .data$badgroup, .data$clinical) %>%
      dplyr::summarise(dplyr::across(columns, ~sum(!is.na(.)))) %>%
      reshape2::melt(id.vars = c("label", "badgroup", "clinical"),
                     measure.vars = columns,
                     variable.name = "classification",
                     value.name = "total")
  } else {
    n <- df %>%
      dplyr::group_by(.data$label, .data$badgroup) %>%
      dplyr::summarise(dplyr::across(columns, ~sum(!is.na(.)))) %>%
      reshape2::melt(id.vars = c("label", "badgroup"),
                     measure.vars = columns,
                     variable.name = "classification",
                     value.name = "total")
  }


  # Count of the number of samples with resistance to each antibiotic class
  # in each hospital-clinical group
  if (by_clinical) {
    s <- df %>%
      dplyr::group_by(.data$label, .data$badgroup, .data$clinical) %>%
      dplyr::summarise(dplyr::across(columns, ~sum(. == 1, na.rm = TRUE))) %>%
      reshape2::melt(id.vars = c("label", "badgroup", "clinical"),
                     measure.vars = columns,
                     variable.name = "classification",
                     value.name = "resistant")
  } else {
    s <- df %>%
      dplyr::group_by(.data$label, .data$badgroup) %>%
      dplyr::summarise(dplyr::across(columns, ~sum(. == 1, na.rm = TRUE))) %>%
      reshape2::melt(id.vars = c("label", "badgroup"),
                     measure.vars = columns,
                     variable.name = "classification",
                     value.name = "resistant")
  }

  # Initialise variables for plotting
  if(this_many <= 12) {
    manual_colours <- stats::setNames(ggthemes::ptol_pal()(this_many),
                                      labs)

  } else {
    colours <- ggthemes::ptol_pal()(12)
    manual_colours <- stats::setNames(
      grDevices::colorRampPalette(colours)(this_many), labs)
  }

  manual_shapes <- stats::setNames(c(16, 16, 1:(this_many-2)), labs)

  labs <- labs[which(labs %in% s$label)]

  # Calculate the proportion of samples in each hospital-clinical group
  # that are resistant to each antibiotic class
  if (by_clinical) {
    probabilities <- merge(n, s, by = c("label", "badgroup", "clinical",
                                        "classification")) %>%
      dplyr::mutate(clinical = dplyr::case_when(
        clinical == "yes" ~ "Clinical",
        clinical == "no" ~ "Carriage"))

  } else {
    probabilities <- merge(n, s, by = c("label", "badgroup", "classification"))
  }

  probabilities <- probabilities %>%
    dplyr::rename(goodbad = .data$badgroup) %>%
    dplyr::mutate(Probability = .data$resistant / .data$total,
                  goodbad = dplyr::case_when(
                    goodbad == 1 ~ "Bad group",
                    goodbad == 0 ~ "Good group")) %>%
    dplyr::mutate(goodbad = factor(.data$goodbad,
                                   levels = c("Good group", "Bad group"))) %>%
    merge(labels$good) %>%
    dplyr::mutate(name = dplyr::case_when(
      grepl("Hospital", .data$label) ~ "Hospital",
      T ~ .data$label),
      label = factor(.data$label),
      label = forcats::fct_relevel(.data$label, labs))
  labs <- levels(probabilities$label)

  # Plot violins

  if (by_clinical) {
    g <- dat %>%
      ggplot2::ggplot() +
      geom_split_violin(ggplot2::aes(x = .data$classification,
                                     y = .data$value,
                                     fill = clinical),
                        trim = TRUE,
                        colour = "grey70",
                        scale = "width")
  } else {
    g <- dat %>%
      ggplot2::ggplot() +
      ggplot2::geom_violin(ggplot2::aes(x = .data$classification,
                                        y = .data$value),
                           trim = TRUE,
                           colour = "grey70",
                           fill = "grey90",
                           scale = "width")
  }

  g <- g + ggplot2::theme_minimal() +
    ggplot2::facet_grid(.~goodbad, scales = "free") +
    ggplot2::scale_fill_manual(name = "",
                               values = c("grey90", "white")) +
    ggplot2::stat_summary(ggplot2::aes(x = .data$classification,
                                       y = .data$value),
                          fun = stats::median,
                          geom = "point",
                          position = ggplot2::position_dodge(width = 1),
                          colour = "grey33") +
    ggplot2::stat_summary(ggplot2::aes(x = .data$classification,
                                       y = .data$value),
                          fun = function(z) {
                            stats::quantile(z, c(0.25, 0.75))
                          },
                          geom = "line",
                          colour = "grey33") +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title = "Plot 2",
                  y = "Probability of resistance") +
    ggplot2::coord_flip(ylim = c(0, 1))

  manual_colours <- manual_colours[labs]
  manual_shapes <- manual_shapes[labs]

  g + ggplot2::geom_point(ggplot2::aes(x = .data$classification,
                                       y = .data$Probability,
                                       colour = .data$label,
                                       shape = .data$label),
                          size = 2,
                          probabilities) +
    ggplot2::scale_colour_manual(name = "",
                                 values = manual_colours,
                                 labels = labs) +
    ggplot2::scale_shape_manual(name = "",
                                values = manual_shapes,
                                labels = labs) +
    ggplot2::theme(legend.position = "bottom")
}


plot3 <- function(model.ggs, parameters) {
  plot_these <- sapply(c("intercept", "sd"), function(x)
    parameters[grepl(x, parameters)]) %>%
    unlist() %>%
    unique() %>%
    c()
  title <- ""

  # Get posteriors
  plot_this <- model.ggs %>%
    dplyr::filter(.data$Parameter %in% plot_these)


  # Generate plot
  g <- plot_this %>%
    ggplot2::ggplot() +
    ggplot2::facet_grid(Parameter~., scale = "free") +
    ggplot2::theme_minimal() +
    ggplot2::geom_violin(ggplot2::aes_string(x = "Parameter",
                                             y = "value"),
                         trim = TRUE,
                         colour = "grey70",
                         fill = "grey90",
                         scale = "width") +
    ggplot2::stat_summary(ggplot2::aes_string(x = "Parameter", y = "value"),
                          fun = stats::median,
                          geom = "point",
                          position = ggplot2::position_dodge(width = 1),
                          colour = "grey33") +
    ggplot2::stat_summary(ggplot2::aes_string(x = "Parameter", y = "value"),
                          fun = function(z) {
                            stats::quantile(z, c(0.25, 0.75))
                          },
                          geom = "line",
                          colour = "grey33") +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title = "Plot 3",
                  y = "value") +
    ggplot2::coord_flip()
}

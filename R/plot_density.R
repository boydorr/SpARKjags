#' plot_density
#'
#' Generates a violin plot where violins represent the posterior probability of
#' resistance and points represent the calculated probability of resistance to
#' each antibiotic class.
#'
#' @param model a \code{runjags} object containing model results
#' @param data data input
#' @param save_to dir
#' @param model_name model name
#' @param save_pdf boolean
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data <- jags_data(classification = "all",
#'                   categories = "human",
#'                   pathogen = "Klebsiella pneumoniae",
#'                   removeQuinPen = TRUE)
#' res.a <- get_model("a", "goodbad_models")
#' plot_density(model = res.a,
#'             data = data)
#' }
#'
plot_density <- function(model,
                         data,
                         save_to,
                         model_name,
                         save_pdf) {
  # Is the plot cached?
  save_to <- file.path(save_to, "density_plots")
  filename <- paste0(gsub("res.", "", model_name), ".rds")
  filepath <- file.path(save_to, filename)

  if(file.exists(filepath)) {
    output <- readRDS(filepath)

  } else {
    # data.frame listing SpARK samples, and their resistances to each antibiotic
    # class, as well as the posterior probability of being in the bad group
    # (mean.p.bad), which defines the badgroup (1 if mean.p.bad > 0.5)
    df <- import_data(model, data)

    # Labels
    if(any(colnames(df) %in% "badgroup")) {
      labels <- get_labels(data)
      params <- get_params()
      var.regex <- get_vars(model)

    }else {
      tmp <- data$lookup$antibiotic_class %>%
        dplyr::mutate(index = paste0("a.prob[", 1:13, "]"))
      labels <- list(tmp, NA, NA)
      params <- list(`probability of resistance` = "a.prob", "intercept", "sd")
      var.regex <- get_vars(model)
    }

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

      # Initialise parameters
      plot_these <- sapply(params[this_set], function(x)
        parameters[grepl(x, parameters)]) %>%
        unlist() %>%
        unique()
      title <- names(params)[[this_set]]
      good <- if_else(grepl("(good group)", title), T, F)
      bad <- if_else(grepl("\\(bad group)", title), T, F)

      # Get posterior probabilities of resistance to each antibiotic class
      plot_this <- model.ggs %>%
        dplyr::filter(.data$Parameter %in% plot_these)

      g <- plot_this %>%
        ggplot2::ggplot() +
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
        ggplot2::labs(title = title)

      if(title == "probability of being in the bad group") {
        labs <- labels[[1]] %>%
          dplyr::filter(.data$index %in% unique(plot_this$Parameter))
        g <- g + ggplot2::scale_x_discrete(breaks = labs$index,
                                           labels = labs$category)
      }


      # Plot prob of being in the badgroup / prob of resistance --------------

      if(grepl("^probability of resistance", title)) {

        # Initialise variables for plotting
        all_labs <- unique(df$label)

        labs <- c("Hospital\n(Carriage)",
                  "Hospital\n(Clinical)",
                  "GP",
                  "Outpatients",
                  "Volunteers")
        labs <- c(labs, all_labs[!all_labs %in% labs])

        this_many <- length(labs)

        # Extract data
        if(bad) df <- df %>% dplyr::filter(.data$badgroup == 1)
        if(good) df <- df %>% dplyr::filter(.data$badgroup == 0)

        # Note that hospital-clinical groups are defined as:
        # Hospital\n(Carriage),
        # Outpatients,
        # Hospital\n(Clinical),
        # Hospital\n(Clinical),
        # Hospital\n(Carriage), and
        # Hospital\n(Clinical)

        # Calculate the observed resistance probability -----------------------

        # Define which columns should be used when the input is a goodbad model
        # or naive (e.g. res.a_naive)
        if(good | bad) {
          columns <- df %>%
            select(-.data$GUID, -.data$name, -.data$hospital, -.data$clinical,
                   -.data$mean.p.bad, -.data$badgroup, -.data$label) %>%
            colnames()
        } else {
          columns <- df %>%
            select(-.data$GUID, -.data$name, -.data$hospital, -.data$clinical,
                   -.data$label) %>%
            colnames()
        }

        # Count of the number of R and S samples in each hospital-clinical group
        # for each antibiotic class (ignore NAs)
        n <- df %>%
          dplyr::group_by(.data$label) %>%
          dplyr::summarise(dplyr::across(columns, ~sum(!is.na(.)))) %>%
          reshape2::melt(id.vars = "label", measure.vars = columns,
                         variable.name = "classification",
                         value.name = "total")

        # Count of the number of samples with resistance to each antibiotic class
        # in each hospital-clinical group
        s <- df %>%
          dplyr::group_by(.data$label) %>%
          dplyr::summarise(dplyr::across(columns, ~sum(. == 1, na.rm = TRUE))) %>%
          reshape2::melt(id.vars = "label", measure.vars = columns,
                         variable.name = "classification",
                         value.name = "resistant")

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
        probabilities <- merge(n, s, by = c("label", "classification")) %>%
          dplyr::mutate(Probability = .data$resistant / .data$total) %>%
          merge(labels[[this_set]]) %>%
          dplyr::mutate(name = dplyr::case_when(
            grepl("Hospital", label) ~ "Hospital",
            T ~ label),
            label = factor(.data$label),
            label = forcats::fct_relevel(.data$label, labs))

        labs <- levels(probabilities$label)

        manual_colours <- manual_colours[labs]
        manual_shapes <- manual_shapes[labs]

        # Extract labels
        xaxis <- labels[[this_set]] %>%
          dplyr::mutate(ind = gsub("a.prob", "", .data$index)) %>%
          dplyr::mutate(classification = dplyr::case_when(
            classification == "Aminoglycoside" ~
              paste("Aminoglycoside", ind),
            classification == "Carbapenem" ~
              paste("Carbapenem", ind),
            classification == "Cephalosporin" ~
              paste("Cephalosporin", ind),
            classification == "Colistin" ~
              paste("Colistin", ind),
            classification == "Fluoroquinolone" ~
              paste("Fluoroquinolone", ind),
            classification == "Fosfomycin" ~
              paste("Fosfomycin", ind),
            classification == "Monobactam" ~
              paste("Monobactam", ind),
            classification == "Nitrofurantoin" ~
              paste("Nitrofurantoin", ind),
            classification == "Penicillin (Penams)" ~
              paste("Penicillin (Penams)", ind),
            classification == "Penicillin Combination" ~
              paste("Penicillin Combination", ind),
            classification == "Tetracycline" ~
              paste("Tetracycline", ind),
            classification == "Trimethoprim" ~
              paste("Trimethoprim", ind),
            classification == "Trimethoprim/Sulfamethoxazole" ~
              paste("Tri/Sul", ind)))
        xaxis.labels <- setNames(xaxis$classification, xaxis$index)

        g <- g + ggplot2::geom_point(
          ggplot2::aes_string(x = "index",
                              y = "Probability",
                              colour = "label",
                              shape = "label"),
          probabilities) +
          ggplot2::scale_colour_manual(name = "",
                                       values = manual_colours,
                                       labels = labs) +
          ggplot2::scale_shape_manual(name = "",
                                      values = manual_shapes,
                                      labels = labs) +
          ggplot2::scale_x_discrete(breaks = plot_these,
                                    labels = xaxis.labels) +
          ggplot2::theme(legend.position = "bottom")
      }

      # Set axis limits
      if(all(grepl("prob", plot_these))) {
        g <- g + ggplot2::coord_flip(ylim = c(0, 1))
      } else {
        g <- g + ggplot2::coord_flip()
      }
    })

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

  return(output)
}

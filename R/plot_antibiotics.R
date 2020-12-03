#' Generate heatmap of resistance to antibiotics
#'
#' Identify bad group samples from bad.p and look at location, AMR profile, and
#' general distribution
#'
#' @param model a \code{runjags} object containing model results
#' @param data a \code{string} specifying the name of the output (an
#' \code{rds} file containing the density plot)
#'
#' @export
#'
plot_antibiotics <- function(model, data) {

  df <- import_data(model, data)
  response <- df %>% dplyr::select(-.data$name, -.data$hospital,
                                   -.data$clinical, -.data$mean.p.bad,
                                   -.data$badgroup, -.data$label)

  # Determine order
  plotthis <- response %>%
    reshape2::melt(id.var = "GUID",
                   variable.name = "class",
                   value.name = "interpretation") %>%
    merge(df %>%
            dplyr::mutate(name = dplyr::case_when(
              .data$name == "Hospital" & .data$clinical == "yes" ~
                "Hospital\nClinical",
              .data$name == "Hospital" & .data$clinical == "no" ~
                "Hospital\nCarriage",
              T ~ name)) %>%
            dplyr::select(.data$GUID, .data$mean.p.bad, .data$name,
                          .data$badgroup)) %>%
    dplyr::filter(.data$badgroup == 1) %>% # bad group
    # Sort samples
    dplyr::group_by(.data$GUID) %>%
    dplyr::mutate(total.guid = sum(.data$interpretation, na.rm = T)) %>%
    dplyr::arrange(.data$name, desc(.data$total.guid)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(GUID = factor(.data$GUID, levels = unique(.data$GUID))) %>%
    # Sort antibiotic classes
    dplyr::group_by(class) %>%
    dplyr::mutate(total.class = sum(.data$interpretation, na.rm = T)) %>%
    dplyr::arrange(desc(.data$total.class)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(class = factor(.data$class, levels = unique(.data$class)))
  classorder <- levels(plotthis$class)

  # Plot results
  plots <- list()
  groups <- c(1, 0)
  plotthis <- lapply(seq_along(groups), function(i) {
    out <- response %>%
      reshape2::melt(id.var = "GUID",
                     variable.name = "class",
                     value.name = "interpretation") %>%
      merge(df %>%
              dplyr::mutate(name = dplyr::case_when(
                .data$name == "Hospital" & .data$clinical == "yes" ~
                  "Hospital Clinical",
                .data$name == "Hospital" & .data$clinical == "no" ~
                  "Hospital Carriage",
                .data$name == "Volunteers" & .data$clinical == "no" ~
                  "Vol",
                T ~ name)) %>%
              dplyr::select(.data$GUID, .data$mean.p.bad, .data$name,
                            .data$badgroup)) %>%
      dplyr::filter(.data$badgroup == groups[i]) %>%
      # Sort samples
      dplyr::group_by(.data$GUID) %>%
      dplyr::mutate(total.guid = sum(.data$interpretation, na.rm = T)) %>%
      dplyr::arrange(.data$name, desc(.data$total.guid)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GUID = factor(.data$GUID, levels = unique(.data$GUID))) %>%
      # Sort antibiotic classes
      dplyr::group_by(class) %>%
      dplyr::mutate(total.class = sum(.data$interpretation, na.rm = T)) %>%
      dplyr::arrange(desc(.data$total.class)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(class = factor(.data$class, levels = classorder)) %>%
      dplyr::group_by(.data$name, .data$badgroup) %>%
      dplyr::mutate(xaxis = as.numeric(.data$GUID)) %>%
      dplyr::mutate(xaxis = .data$xaxis - min(.data$xaxis)) %>%
      dplyr::ungroup()
  })
  plotthis <- do.call(rbind.data.frame, plotthis) %>%
    dplyr::mutate(interpretation = dplyr::case_when(
      .data$interpretation == 1 ~ "Resistant",
      .data$interpretation == 0 ~ "Susceptible"),
      badgroup = dplyr::case_when(.data$badgroup == 1 ~ "Bad group",
                                  .data$badgroup == 0 ~ "Good group"))

  # Plot
  breaks <- seq(0, 500, by = 10) - 0.5
  labels <- c(0, rep("", 4), 50, rep("", 4), 100, rep("", 4), 150, rep("", 4),
              200, rep("", 4), 250, rep("", 4), 300, rep("", 4), 350,
              rep("", 4), 400, rep("", 4), 450, rep("", 4), 500)
  output <- plotthis %>% ggplot2::ggplot() + ggplot2::theme_minimal() +
    ggplot2::geom_tile(ggplot2::aes_string(x = "class", y = "xaxis",
                                           fill = "interpretation"),
                       colour = "grey") + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("#d73027", "#fee08b")) +
    ggplot2::facet_grid(badgroup ~ name, scales = "free", space = "free") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                breaks = breaks,
                                labels = labels) +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_rect(size = 0.5, colour = "grey",
                                               fill = NA),
                   legend.position = "bottom",
                   axis.ticks.x = element_line(size = 0.25, colour = "grey")) +
    ggplot2::labs(x = "Antibiotic class", y = "SpARK sample", fill = "")

  output
}

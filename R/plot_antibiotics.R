#' plot_antibiotics
#'
#' Identify bad group samples from bad.p and look at location, amr profile, and
#' general distribution
#'
#' @param model model
#' @param data data
#'
#' @export
#'
plot_antibiotics <- function(model, data) {

  df <- import_data(model, data)
  response <- df %>% dplyr::select(-.data$name, -.data$hospital,
                                   -.data$clinical, -.data$mean.p.bad,
                                   -.data$badgroup, -.data$label)

  # Plot results
  plots <- list()
  groups <- c(1, 0)
  for (i in seq_along(groups)) {
    plotthis <- response %>%
      dplyr::rename(Ami = .data$Aminoglycoside,
                    Carbape = .data$Carbapenem,
                    Cephalo = .data$Cephalosporin,
                    Colisti = .data$Colistin,
                    Fluoroq = .data$Fluoroquinolone,
                    Fosfomy = .data$Fosfomycin,
                    Monobac = .data$Monobactam,
                    Nitrofu = .data$Nitrofurantoin,
                    `Pen (P)` = .data$`Penicillin (Penams)`,
                    `Pen (C)` = .data$`Penicillin Combination`,
                    Tetracy = .data$Tetracycline,
                    Trimeth = .data$Trimethoprim,
                    `Tri/Sul` = .data$`Trimethoprim/Sulfamethoxazole`) %>%
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
      dplyr::filter(.data$badgroup == groups[i]) %>%
      # Sort samples
      dplyr::group_by(.data$GUID) %>%
      dplyr::mutate(total.guid = sum(.data$interpretation, na.rm = T)) %>%
      dplyr::arrange(.data$name, .data$total.guid) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GUID = factor(.data$GUID, levels = unique(.data$GUID))) %>%
      # Sort antibiotic classes
      dplyr::group_by(class) %>%
      dplyr::mutate(total.class = sum(.data$interpretation, na.rm = T)) %>%
      dplyr::arrange(desc(.data$total.class)) %>%
      dplyr::ungroup()

    if(groups[i] == 1) {
      plotthis <- plotthis %>%
        dplyr::mutate(class = factor(.data$class, levels = unique(.data$class)))
      classorder <- levels(plotthis$class)
      title <- "Bad group"
    }

    if(groups[i] == 0) {
      plotthis <- plotthis %>%
        dplyr::mutate(class = factor(.data$class, levels = classorder))
      title <- "Good group"
    }

    # Tidy labels
    plotthis <- plotthis %>%
      dplyr::mutate(interpretation = dplyr::case_when(
        .data$interpretation == 1 ~ "Resistant",
        .data$interpretation == 0 ~ "Susceptible"),
        badgroup = dplyr::case_when(.data$badgroup == 1 ~ "Bad group",
                                    .data$badgroup == 0 ~ "Good group"))

    # Plot
    g <- plotthis %>% ggplot2::ggplot() + ggplot2::theme_minimal() +
      ggplot2::geom_tile(ggplot2::aes_string(x = "class", y = "GUID",
                                             fill = "interpretation"),
                         colour = "grey") +
      ggplot2::scale_fill_manual(values = c("#d73027", "#fee08b")) +
      ggplot2::facet_grid(name ~ ., scales = "free", space = "free") +
      ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
      ggplot2::labs(x = "Antibiotic class", y = "SpARK ID", title = title) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     strip.text.y = ggplot2::element_text(angle = 0)) +
      ggplot2::labs(fill = "") +
      ggplot2::theme(legend.position = "bottom")

    plots[[i]] <- g
  }

  egg::ggarrange(plots = plots, ncol = 2)
}

#' antibioticsplot
#'
#' Identify bad group samples from bad.p and look at location, amr profile, and
#' general distribution
#' @export
#'
antibioticsplot <- function(model, data) {

  df <- import_data(model, data)
  response <- df %>% dplyr::select(-name, -hospital, -clinical,
                                   -mean.p.bad, -badgroup, -label)

  # Plot results
  plots <- list()
  groups <- c(1, 0)
  for (i in seq_along(groups)) {
    plotthis <- response %>%
      dplyr::rename(Ami = Aminoglycoside,
                    Carbape = Carbapenem,
                    Cephalo = Cephalosporin,
                    Colisti = Colistin,
                    Fluoroq = Fluoroquinolone,
                    Fosfomy = Fosfomycin,
                    Monobac = Monobactam,
                    Nitrofu = Nitrofurantoin,
                    `Pen (P)` = `Penicillin (Penams)`,
                    `Pen (C)` = `Penicillin Combination`,
                    Tetracy = Tetracycline,
                    Trimeth = Trimethoprim,
                    `Tri/Sul` = `Trimethoprim/Sulfamethoxazole`) %>%
      reshape2::melt(id.var = "GUID",
                     variable.name = "class",
                     value.name = "interpretation") %>%
      merge(df %>%
              dplyr::mutate(name = dplyr::case_when(
                name == "Hospital" & clinical == "yes" ~ "Hospital\nClinical",
                name == "Hospital" & clinical == "no" ~ "Hospital\nCarriage",
                T ~ name)) %>%
              dplyr::select(GUID, mean.p.bad, name, badgroup)) %>%
      dplyr::filter(badgroup == groups[i]) %>%
      # Sort samples
      dplyr::group_by(GUID) %>%
      dplyr::mutate(total.guid = sum(interpretation, na.rm = T)) %>%
      dplyr::arrange(name, total.guid) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GUID = factor(GUID, levels = unique(GUID))) %>%
      # Sort antibiotic classes
      dplyr::group_by(class) %>%
      dplyr::mutate(total.class = sum(interpretation, na.rm = T)) %>%
      dplyr::arrange(desc(total.class)) %>%
      dplyr::ungroup()

    if(groups[i] == 1) {
      plotthis %<>%
        dplyr::mutate(class = factor(class, levels = unique(class)))
      classorder <- levels(plotthis$class)
      title <- "Bad group"
    }

    if(groups[i] == 0) {
      plotthis %<>%
        dplyr::mutate(class = factor(class, levels = classorder))
      title <- "Good group"
    }

    # Tidy labels
    plotthis %<>%
      dplyr::mutate(interpretation = dplyr::case_when(
        interpretation == 1 ~ "Resistant",
        interpretation == 0 ~ "Susceptible"),
        badgroup = dplyr::case_when(badgroup == 1 ~ "Bad group",
                                    badgroup == 0 ~ "Good group"))

    # Plot
    g <- plotthis %>% ggplot2::ggplot() + ggplot2::theme_minimal() +
      ggplot2::geom_tile(ggplot2::aes(x = class, y = GUID,
                                      fill = interpretation), colour = "grey") +
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

  cowplot::plot_grid(plotlist = plots, ncol = 2)
}



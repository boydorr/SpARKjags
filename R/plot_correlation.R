#' plot_correlation
#'
#' @param model a \code{runjags} object containing model results
#' @param data a \code{list} containing the data input for the runjags model,
#' generated with \code{jags_data()}
#'
#' @export
#'
plot_correlation <- function(model, data) {

  df <- import_data(model, data)
  response <- df %>%
    dplyr::select(-.data$name, -.data$hospital, -.data$clinical,
                  -.data$mean.p.bad, -.data$badgroup, -.data$label)

  # Plot results
  plots <- list()
  groups <- c(1, 0)
  for (i in seq_along(groups)) {
    these <- response %>%
      reshape2::melt(id.var = "GUID",
                     variable.name = "class",
                     value.name = "interpretation") %>%
      merge(df %>%
              dplyr::mutate(name = dplyr::case_when(
                name == "Hospital" & clinical == "yes" ~ "Hospital\nClinical",
                name == "Hospital" & clinical == "no" ~ "Hospital\nCarriage",
                T ~ name)) %>%
              dplyr::select(.data$GUID, .data$mean.p.bad, .data$name,
                            .data$badgroup)) %>%
      dplyr::filter(.data$badgroup == groups[i]) %>%
      dplyr::select(.data$GUID, .data$badgroup) %>%
      unique()
    these <- these$GUID

    tmp <- response %>%
      dplyr::filter(.data$GUID %in% these) %>%
      dplyr::select(-.data$GUID) %>%
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
                    `Tri/Sul` = .data$`Trimethoprim/Sulfamethoxazole`)

    res <- stats::cor(tmp,
               use = "pairwise.complete.obs",
               method = "pearson")

    g <- GGally::ggcorr(res,
                        hjust = 0.75,
                        label = T,
                        label_alpha = T)

    if(groups[i] == 1) g <- g + ggplot2::labs(title = "Bad group")
    if(groups[i] == 0) g <- g + ggplot2::labs(title = "Good group")

    plots[[i]] <- g
  }

  egg::ggarrange(plots = plots, nrow = 2)
}


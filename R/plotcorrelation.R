#' plotcorrelation
#'
#' @export
#'
plotcorrelation <- function(model, data) {

  df <- import_data(model, data)
  response <- df %>% dplyr::select(-name, -hospital, -clinical,
                                   -mean.p.bad, -badgroup, -label)

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
              dplyr::select(GUID, mean.p.bad, name, badgroup)) %>%
      dplyr::filter(badgroup == groups[i]) %>%
      dplyr::select(GUID, badgroup) %>%
      unique() %$%
      GUID

    tmp <- response %>%
      dplyr::filter(GUID %in% these) %>%
      dplyr::select(-GUID) %>%
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
                    `Tri/Sul` = `Trimethoprim/Sulfamethoxazole`)

    res <- cor(tmp,
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

  cowplot::plot_grid(plotlist = plots, ncol = 1)
}







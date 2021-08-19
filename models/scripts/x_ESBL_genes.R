#' ---
#' title: "ESBL gene summary"
#' author: "Sonia Mitchell"
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#' ---


#+ r setup, include=FALSE, echo = F

library(SpARK)
library(dplyr)


#+

#' # ESBL genes {.tabset}
#' ## KLEBdata column names
colnames(KLEBdata)

#' #' ## Unique genes
#' apply(KLEBdata, 2, function(x) length(unique(x)))
#' which(colnames(KLEBdata) %in% "AGly")
#' apply(KLEBdata[80:99], 2, unique)
#'

#' ## Unique ESBL genes
data <- KLEBdata %>%
  dplyr::filter(species == "Klebsiella pneumoniae") %>%
  dplyr::select(GUID, Bla_ESBL) %>%
  dplyr::filter(Bla_ESBL != "-") %>%
  dplyr::mutate(Bla_ESBL = dplyr::if_else(grepl(";", Bla_ESBL),
                                          Bla_ESBL, paste0(Bla_ESBL, ";"))) %>%
  tidyr::separate(Bla_ESBL, c("gene1", "gene2"), sep = ";")

data_long <- dplyr::bind_rows(data %>%  dplyr::select(GUID, gene1),
                              data %>%  dplyr::select(GUID, gene2) %>%
                                dplyr::rename(gene1 = gene2)) %>%
  dplyr::filter(gene1 != "") %>%
  dplyr::rename(gene = gene1)

data_long$gene %>% unique() %>% sort()

#' ## ESBL gene plot
meta <- METAdata %>%
  dplyr::filter(GUID %in% data_long$GUID) %>%
  dplyr::select(GUID, ASSOCIATED_SPECIES, ASSOCIATED_GROUP, TYPE, Category) %>%
  dplyr::mutate(tag = dplyr::case_when(Category == "human" ~ TYPE,
                                       Category == "other" ~ TYPE,
                                       T ~ ASSOCIATED_SPECIES))

plot_this <- merge(data_long, meta) %>%
  dplyr::group_by(gene, tag) %>%
  dplyr::summarise(total = n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(gene = unique(.$gene),
                  tag = unique(.$tag),
                  fill = list(total = NA)) %>%
  dplyr::arrange(gene) %>%
  dplyr::mutate(gene = factor(gene, levels = unique(gene)))

total_count <- plot_this %>%
  dplyr::group_by(tag) %>%
  dplyr::summarise(total_count = sum(total, na.rm = T))

labs <- plot_this %>%
  dplyr::select(tag) %>%
  unique() %>%
  dplyr::mutate(tag_numeric = as.numeric(as.factor(tag))) %>%
  merge(total_count)

plot_this %<>% merge(labs)




plot_this %>% ggplot2::ggplot() + ggplot2::theme_minimal() +
  ggplot2::geom_tile(ggplot2::aes(x = tag_numeric, y = gene, fill = total),
                     colour = "grey") +
  ggplot2::scale_fill_gradient(name = "Count",
                               low = "#fff7bc",
                               high = "#d95f0e",
                               na.value = "#a1d99b") +
  ggplot2::geom_text(ggplot2::aes(x = tag_numeric, y = gene, label = total),
                     color = "grey30") +
  ggplot2::labs(x = "Number of ESBL gene occurances",
                y = "ESBL genes") +
  # ggplot2::scale_x_discrete(position = "top") +
  ggplot2::scale_x_continuous(breaks = labs$tag_numeric,
                              labels = labs$total_count,
                              sec.axis = ggplot2::sec_axis(
                                ~.,
                                breaks = labs$tag_numeric,
                                labels = labs$tag)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0))














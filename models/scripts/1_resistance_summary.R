#' ---
#' title: Jags models
#' date: Sep 2019
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#' ---
#' Last updated: `r Sys.time()`

#+ r setup, include=FALSE

knitr::opts_chunk$set(message=FALSE)
library(SpARKcarbapenam)
library(SpARK)
library(dplyr)

save_plots <- FALSE

#+

df <- METAdata %>%
  merge(ATBdata %>%
          dplyr::rename(GUID = UNIQUE_SPARK_ID)) %>%
  dplyr::filter(used_MIC == "yes",
                !grepl("NEG", GUID),
                DATE_OF_BIRTH != "XXXX") %>%
  dplyr::select(GUID, Category, GROUP, TYPE, Clinical, SPECIFIC_GROUP,
                SAMPLE_TYPE, Interpretation, Antibiotic_name, Classification,
                Phoenix_Organism) %>%
  dplyr::filter(Interpretation == "R")

df_kleb <- df %>%
  merge(KLEBdata %>%
          dplyr::select(GUID, virulence_score, resistance_score))

labs <- cbind.data.frame(Classification = unique(df$Classification),
                         tag = substr(unique(df$Classification), 1, 3),
                         stringsAsFactors = F) %>%
  dplyr::mutate(tag = dplyr::case_when(
    Classification == "Aminoglycoside" ~ "Aminoglycoside",
    Classification == "Cephalosporin" ~ "Cephalosporin",
    Classification == "Carbapenem" ~ "Carbapenem",
    Classification == "Fluoroquinolone" ~ "Fluoroquinolone",
    Classification == "Penicillin (Penams)" ~ "Pen P",
    Classification == "Trimethoprim/Sulfamethoxazole" ~ "Tri/Sul",
    Classification == "Penicillin Combination" ~ "Pen C",
    Classification == "Cephalosporin" ~ "Ceph",
    T ~ tag))


# All samples -------------------------------------------------------------

# All antibiotic classes; human category
all_human <- df %>%
  filter(Category == "human") %>%
  merge(labs) %>%
  mutate(Clinical = dplyr::case_when(
    SPECIFIC_GROUP == "Pavia_citizen" ~ "no",
    T ~ Clinical),
    facet = dplyr::case_when(
      SPECIFIC_GROUP == "Pavia_citizen" ~ "Volunteer\n(Carriage)",
      SPECIFIC_GROUP == "Belgioioso" ~ "GP\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "yes" ~
        "Hospital\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "no" ~
        "Hospital\n(Carriage)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "dontknow" ~
        "Hospital\n(Don't know)"),
    facet = factor(facet, levels = c("Hospital\n(Clinical)",
                                     "Hospital\n(Carriage)",
                                     "Hospital\n(Don't know)",
                                     "GP\n(Clinical)",
                                     "Volunteer\n(Carriage)")),
    col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                    T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes; animal category
all_animal <- df %>%
  filter(Category == "animal") %>%
  merge(labs) %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes; livestock
all_livestock <- df %>%
  filter(Category == "animal",
         TYPE == "livestock") %>%
  merge(labs) %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes; other category
all_other <- df %>%
  filter(Category == "other") %>%
  merge(labs) %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

#' ## All antibiotic classes

#+ fig.height = 10, fig.width = 8
all_data <- rbind.data.frame(all_human, all_livestock,
                             all_animal, all_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(tag~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = Antibiotic_name, y = Count, fill = col),
                    stat = "identity", colour = "black") +
  ggplot2::geom_text(ggplot2::aes(x = Antibiotic_name, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_fill_manual(name = "", values = c("deeppink4", "gray66")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .3))) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Antibiotic", y = "Number of resistant samples")
all_data

if(save_plots)
  ggplot2::ggsave('alldata.pdf', all_data, width = 14, height = 12)





# Remove samples with carbapenem resistance -------------------------------

# Samples with carbapenem resistance
samples <- df %>%
  dplyr::filter(Classification == "Carbapenem" ,
                Interpretation == "R") %$%
  GUID %>%
  unique()

df_nocarb <- df %>%
  dplyr::filter(!GUID %in% samples)

# All antibiotic classes except carbapenem; human category
nc_human <- df_nocarb %>%
  filter(Category == "human") %>%
  merge(labs) %>%
  dplyr::mutate(Clinical = dplyr::case_when(
    SPECIFIC_GROUP == "Pavia_citizen" ~ "no",
    T ~ Clinical),
    facet = dplyr::case_when(
      SPECIFIC_GROUP == "Pavia_citizen" ~ "Volunteer\n(Carriage)",
      SPECIFIC_GROUP == "Belgioioso" ~ "GP\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "yes" ~
        "Hospital\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "no" ~
        "Hospital\n(Carriage)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "dontknow" ~
        "Hospital\n(Don't know)"),
    facet = factor(facet, levels = c("Hospital\n(Clinical)",
                                     "Hospital\n(Carriage)",
                                     "Hospital\n(Don't know)",
                                     "GP\n(Clinical)",
                                     "Volunteer\n(Carriage)")),
    col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                    T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes except carbapenem; animal category
nc_animal <- df_nocarb %>%
  filter(Category == "animal",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes except carbapenem; livestock
nc_livestock <- df_nocarb %>%
  filter(Category == "animal",
         TYPE == "livestock",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# All antibiotic classes except carbapenem; other category
nc_other <- df_nocarb %>%
  filter(Category == "other",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

#' ## All antibiotic classes except carbapenem

#+ fig.height = 10, fig.width = 8
nc_data <- rbind.data.frame(nc_human, nc_livestock,
                            nc_animal, nc_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(tag~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = Antibiotic_name, y = Count, fill = col),
                    stat = "identity", colour = "black", fill = "gray66") +
  ggplot2::geom_text(ggplot2::aes(x = Antibiotic_name, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .3))) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Antibiotic", y = "Number of resistant samples")
nc_data

if(save_plots)
  ggplot2::ggsave('ncdata.pdf', nc_data, width = 14, height = 12)


# Remove samples without carbapenem resistance ----------------------------

# Samples with carbapenem resistance
samples <- df %>%
  dplyr::filter(Classification == "Carbapenem" ,
                Interpretation == "R") %$%
  GUID %>%
  unique()

df_carb <- df %>%
  dplyr::filter(GUID %in% samples)

# Samples with carbapenem resistance; human category
nc_human <- df_carb %>%
  filter(Category == "human") %>%
  merge(labs) %>%
  dplyr::mutate(Clinical = dplyr::case_when(
    SPECIFIC_GROUP == "Pavia_citizen" ~ "no",
    T ~ Clinical),
    facet = dplyr::case_when(
      SPECIFIC_GROUP == "Pavia_citizen" ~ "Volunteer\n(Carriage)",
      SPECIFIC_GROUP == "Belgioioso" ~ "GP\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "yes" ~
        "Hospital\n(Clinical)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "no" ~
        "Hospital\n(Carriage)",
      SPECIFIC_GROUP %in% c("Montescano", "Maugeri", "Santa_Margherita",
                            "San_Matteo") & Clinical == "dontknow" ~
        "Hospital\n(Don't know)"),
    facet = factor(facet, levels = c("Hospital\n(Clinical)",
                                     "Hospital\n(Carriage)",
                                     "Hospital\n(Don't know)",
                                     "GP\n(Clinical)",
                                     "Volunteer\n(Carriage)")),
    col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                    T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# Samples with carbapenem resistance; animal category
nc_animal <- df_carb %>%
  filter(Category == "animal",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# Samples with carbapenem resistance; livestock
nc_livestock <- df_carb %>%
  filter(Category == "animal",
         TYPE == "livestock",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

# Samples with carbapenem resistance; other category
nc_other <- df_carb %>%
  filter(Category == "other",
         Classification != "Carbapenem") %>%
  merge(labs) %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(tag, Antibiotic_name, facet, col) %>%
  summarise(Count = n())

#' ## Samples with carbapenem resistance

#+ fig.height = 10, fig.width = 8
c_data <- rbind.data.frame(nc_human, nc_livestock,
                           nc_animal, nc_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(tag~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = Antibiotic_name, y = Count, fill = col),
                    stat = "identity", colour = "black") +
  ggplot2::geom_text(ggplot2::aes(x = Antibiotic_name, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .3))) +
  ggplot2::scale_fill_manual(name = "", values = c("deeppink4", "gray66")) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Antibiotic", y = "Number of resistant samples")
c_data

if(save_plots)
  ggplot2::ggsave('cdata.pdf', c_data, width = 14, height = 12)






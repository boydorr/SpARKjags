# generate all.data style plots with all data, carbapenem samples, samples without carbapenem resistance, but instead of antibiotic class, put virulence class on the y axis

# filter out carbapenem samples
# run jags models on other antibiotics
# do models look similar?
# is clinical / carriage still important?


library(SpARKcarbapenam)
library(SpARK)
library(dplyr)
# library(tidyr)
# library(ggplot2)

#+

unique(METAdata$ASSOCIATED_SPECIES) %>% sort()


df <- METAdata %>%
  merge(ATBdata %>%
          dplyr::rename(GUID = UNIQUE_SPARK_ID)) %>%
  dplyr::filter(used_MIC == "yes",
                !grepl("NEG", GUID),
                DATE_OF_BIRTH != "XXXX") %>%
  dplyr::select(GUID, Category, GROUP, TYPE, Clinical, SPECIFIC_GROUP,
                SAMPLE_TYPE, Interpretation, Antibiotic_name, Classification,
                Phoenix_Organism, ASSOCIATED_SPECIES)

df_kleb <- df %>%
  merge(KLEBdata %>%
          dplyr::select(GUID, virulence_score, resistance_score))



# All samples -------------------------------------------------------------

# All antibiotic classes; human category
all_human <- df_kleb %>%
  filter(Category == "human") %>%
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
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; animal category
all_animal <- df_kleb %>%
  filter(Category == "animal") %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

df_kleb %>%
  filter(Category == "animal") %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet, ASSOCIATED_SPECIES) %>%
  summarise(Count = n_distinct(GUID)) %>%
  arrange(desc(Count), desc(virulence_score))


# All antibiotic classes; livestock
all_livestock <- df_kleb %>%
  filter(Category == "animal",
         TYPE == "livestock") %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cattle
all_cattle <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cattle") %>%
  mutate(facet = "cattle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; pigs
all_pigs <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "pig") %>%
  mutate(facet = "pig",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; dogs
all_dogs <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "dog") %>%
  mutate(facet = "dog",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cats
all_cats <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cat") %>%
  mutate(facet = "cat",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; turtles
all_turtles <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "turtle") %>%
  mutate(facet = "turtle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; flies
all_flies <- df_kleb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "fly") %>%
  mutate(facet = "fly",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; other category
all_other <- df_kleb %>%
  filter(Category == "other") %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

all_data <- rbind.data.frame(all_human, all_livestock, all_animal,
                             all_cattle, all_pigs, all_dogs, all_cats,
                             all_turtles, all_flies, all_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(.~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = virulence_score, y = Count),
                    stat = "identity", colour = "black", fill = "gray66") +
  ggplot2::geom_text(ggplot2::aes(x = virulence_score, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(0, .4))) +
  ggplot2::scale_x_continuous(breaks = 0:5) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Virulence score", y = "Number of resistant samples",
                title = "All samples")
ggplot2::ggsave('vir_alldata.pdf', all_data, width = 24, height = 4)





# Remove samples with carbapenem resistance -------------------------------

# Samples with carbapenem resistance
samples <- df_kleb %>%
  dplyr::filter(Classification == "Carbapenem" ,
                Interpretation == "R") %$%
  GUID %>%
  unique()

df_nocarb <- df_kleb %>%
  dplyr::filter(!GUID %in% samples)

# All antibiotic classes except carbapenem; human category
nc_human <- df_nocarb %>%
  filter(Category == "human") %>%
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
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; animal category
nc_animal <- df_nocarb %>%
  filter(Category == "animal",
         Classification != "Carbapenem") %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; livestock
nc_livestock <- df_nocarb %>%
  filter(Category == "animal",
         TYPE == "livestock",
         Classification != "Carbapenem") %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cattle
nc_cattle <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cattle") %>%
  mutate(facet = "cattle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; pigs
nc_pigs <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "pig") %>%
  mutate(facet = "pig",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; dogs
nc_dogs <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "dog") %>%
  mutate(facet = "dog",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cats
nc_cats <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cat") %>%
  mutate(facet = "cat",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; turtles
nc_turtles <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "turtle") %>%
  mutate(facet = "turtle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; flies
nc_flies <- df_nocarb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "fly") %>%
  mutate(facet = "fly",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; other category
nc_other <- df_nocarb %>%
  filter(Category == "other",
         Classification != "Carbapenem") %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

nc_data <- rbind.data.frame(nc_human, nc_livestock, nc_animal,
                            nc_cattle, nc_pigs, nc_dogs, nc_cats,
                            nc_turtles, nc_flies, nc_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(.~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = virulence_score, y = Count, fill = col),
                    stat = "identity", colour = "black", fill = "gray66") +
  ggplot2::geom_text(ggplot2::aes(x = virulence_score, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(0, .4))) +
  ggplot2::scale_x_continuous(breaks = 0:5) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Virulence score", y = "Number of resistant samples",
                title = "No carbapenam resistance")

ggplot2::ggsave('vir_ncdata.pdf', nc_data, width = 24, height = 4)


# Remove samples without carbapenem resistance ----------------------------

# Samples with carbapenem resistance
samples <- df_kleb %>%
  dplyr::filter(Classification == "Carbapenem" ,
                Interpretation == "R") %$%
  GUID %>%
  unique()

df_carb <- df_kleb %>%
  dplyr::filter(GUID %in% samples)

# All antibiotic classes except carbapenem; human category
nc_human <- df_carb %>%
  filter(Category == "human") %>%
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
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; animal category
nc_animal <- df_carb %>%
  filter(Category == "animal",
         Classification != "Carbapenem") %>%
  mutate(facet = "animal",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; livestock
nc_livestock <- df_carb %>%
  filter(Category == "animal",
         TYPE == "livestock",
         Classification != "Carbapenem") %>%
  mutate(facet = "livestock",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cattle
nc_cattle <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cattle") %>%
  mutate(facet = "cattle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; pigs
nc_pigs <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "pig") %>%
  mutate(facet = "pig",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; dogs
nc_dogs <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "dog") %>%
  mutate(facet = "dog",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; cats
nc_cats <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "cat") %>%
  mutate(facet = "cat",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; turtles
nc_turtles <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "turtle") %>%
  mutate(facet = "turtle",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; flies
nc_flies <- df_carb %>%
  filter(Category == "animal",
         ASSOCIATED_SPECIES == "fly") %>%
  mutate(facet = "fly",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

# All antibiotic classes; other category
nc_other <- df_carb %>%
  filter(Category == "other",
         Classification != "Carbapenem") %>%
  mutate(facet = "other",
         col = case_when(Classification == "Carbapenem" ~ "Carbapenem",
                         T ~ "Other")) %>%
  group_by(virulence_score, facet) %>%
  summarise(Count = n_distinct(GUID))

c_data <- rbind.data.frame(nc_human, nc_livestock, nc_animal,
                           nc_cattle, nc_pigs, nc_dogs, nc_cats,
                           nc_turtles, nc_flies, nc_other) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::facet_grid(.~facet, scales = "free", space = "free_y") +
  ggplot2::geom_bar(ggplot2::aes(x = virulence_score, y = Count),
                    stat = "identity", colour = "black", fill = "gray66") +
  ggplot2::geom_text(ggplot2::aes(x = virulence_score, y = Count, label = Count),
                     position = ggplot2::position_dodge(width = .9), hjust = -.2) +
  ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(0, .4))) +
  ggplot2::scale_x_continuous(breaks = 0:5) +
  ggplot2::theme(panel.spacing.y = ggplot2::unit(0,"lines"),
                 panel.border = ggplot2::element_rect(color = "grey60"),
                 strip.placement = "outside",
                 legend.position = "none") +
  ggplot2::labs(x = "Virulence score", y = "Number of resistant samples",
                title = "Samples with carbapenam resistance")

ggplot2::ggsave('vir_cdata.pdf', c_data, width = 24, height = 4)






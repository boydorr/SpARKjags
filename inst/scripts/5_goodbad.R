#' ---
#' title: Goodbad models
#' date: Feb 2020
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#' ---
#' Last updated: `r Sys.time()`

#+ setup, include=F

library(SpARKjags)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE)

directory <- "goodbad_models"

data_full <- jags_data(classification = "all",
                       categories = "human",
                       pathogen = "Klebsiella pneumoniae",
                       removeQuinPen = F)

# Remove Quinolone and Penicillin as levels in antibiotic_classes
data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

data_ah <- jags_data(classification = "all",
                     categories = c("human", "animal"),
                     pathogen = "Klebsiella pneumoniae",
                     removeQuinPen = T)


#+


#' ### Data summary
data_full$jags$response %>%
  data.frame() %>%
  dplyr::mutate_if(is.numeric, as.factor) %>%
  summary() %>%
  data.frame() %>%
  dplyr::select(-Var1) %>%
  tidyr::separate(Freq, c("level", "count"), sep = ":") %>%
  dplyr::mutate(count = as.numeric(count),
                level = gsub(" ", "", level),
                fill = if_else(Var2 %in% c("Quinolone", "Penicillin"),
                               T, F)) %>%
  dplyr::filter(!is.na(level)) %>%
  ggplot2::ggplot() + ggplot2::theme_minimal() +
  ggplot2::geom_bar(ggplot2::aes(x = level, y = count, fill = fill),
                    stat = "identity", colour = "black") +
  ggplot2::scale_fill_manual(values = c("grey", "red")) +
  ggplot2::facet_wrap(~Var2) +
  ggplot2::theme(legend.position = "none")

#' Remove Quinolone and Penicillin as levels in antibiotic_classes
data$data.human$data %>%
  dplyr::group_by(dataset, clinical) %>%
  dplyr::summarise(count = n()) %>%
  data.frame() %>%
  flextable::regulartable() %>%
  flextable::autofit()



#'
#' # 1. Control for the structure of the data
#'
#' * Posterior: The grey bar represents the 1st and 3rd quantile distance around
#'   the median
#' * AMR Correlation: calculated Pearson correlation using pairwise.complete.obs
#' * Class interpretation: R when sample is resistant to at least one
#'   antibiotic, S when sample is susceptible to at least one antibiotic,
#'   otherwise NA
#' * Samples are determined as being in the bad group when bad.p (bad.gp, bad,v,
#'   and bad.o) is greater than 0.5
#'
#' ## 1a. Experimental structure
#'
#' ### res.a_naive {.tabset}
#' response ~ antibiotic_class (naive model)
#'

path <- run_SpARKjags_model(data, file.path(directory, "a_naive.R"))
res.a_naive <- get_model(path)

#' #### Posterior
#+ res.a_naive, fig.height = 6
res.a_naive %>% plot_density(data)

#' #### Diagnostics
res.a_naive %>% DIC() # 6166.735
res.a_naive %>% testSSEF()
res.a_naive %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_naive %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_naive %>% plot_autocorr()

#' #### Model
res.a_naive$model

#' #### Results
res.a_naive




#' ### res.a {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data, file.path(directory, "a.R"), thin = 10)
res.a <- get_model(path)

#' #### Posterior
#+ res.a, fig.height = 10
res.a %>% plot_density(data)

#' #### Statistics
res.a %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.a %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.a %>% plot_correlation(data)

#' #### Diagnostics
res.a %>% DIC() # 6151.664
res.a %>% testSSEF()
res.a %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a %>% plot_autocorr()

#' #### Model
res.a$model

#' #### Results
res.a

data_ah$data.animal$response %>% data.frame() %>%
  select(Carbapenem) %>% sum(na.rm = T)


#' ### Animals {-}


# Full data
SpARK::METAdata %>%
  dplyr::select(GUID, ASSOCIATED_SPECIES, TYPE) %>%
  merge(SpARK::KLEBdata %>% dplyr::select(GUID, species)) %>%
  dplyr::filter(!is.na(ASSOCIATED_SPECIES),
                !grepl("^Un", ASSOCIATED_SPECIES),
                TYPE != "Failed",
                TYPE != "environment",
                ASSOCIATED_SPECIES != "human",
                species == "Klebsiella pneumoniae") %>%
  dplyr::group_by(ASSOCIATED_SPECIES, species, TYPE) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::coord_flip() +
  ggplot2::geom_bar(ggplot2::aes(x = ASSOCIATED_SPECIES, y = count,
                                 fill = TYPE),
                    colour = "black",
                    stat = "identity") +
  ggplot2::facet_wrap(~species, scales = "free")


# Jags data
livestock <- data_ah$data.animal$livestock_dat %>%
  dplyr::group_by(associated_species) %>%
  dplyr::rename(index = associated_species) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  merge(data_ah$data.animal$lookup_tables$associated_species) %>%
  dplyr::mutate(type = "Livestock") %>%
  dplyr::select(-index)

companion <- data_ah$data.animal$companion_dat %>%
  dplyr::group_by(associated_species) %>%
  dplyr::rename(index = associated_species) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  merge(data_ah$data.animal$lookup_tables$associated_species) %>%
  dplyr::mutate(type = "Companion animals") %>%
  dplyr::select(-index)

wild <- data_ah$data.animal$wild_dat %>%
  dplyr::group_by(associated_species) %>%
  dplyr::rename(index = associated_species) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  merge(data_ah$data.animal$lookup_tables$associated_species) %>%
  dplyr::mutate(type = "Wild animals") %>%
  dplyr::select(-index)

dplyr::bind_rows(livestock, companion, wild) %>%
  dplyr::arrange(count) %>%
  dplyr::mutate(associated_species = factor(associated_species,
                                            levels = associated_species)) %>%
  ggplot2::ggplot() + ggplot2::theme_minimal() + ggplot2::coord_flip() +
  ggplot2::geom_bar(ggplot2::aes(x = associated_species, y = count,
                                 fill = type),
                    colour = "black", stat = "identity") +
  ggplot2::geom_text(ggplot2::aes(x = associated_species, y = count,
                                  label = count),
                     position = ggplot2::position_dodge(width = 0.9),
                     hjust = -0.25) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .1))) +
  ggplot2::labs(title = "Animal samples", fill = "", x = "Associated species",
                y = "Number of samples") +
  ggplot2::scale_fill_manual(values = c("#33658A", "#2F4858", "#F6AE2D")) +
  ggplot2::theme(legend.position = c(0.8, 0.5))




#' ### res.a_livestock {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_livestock.R"),
                            thin = 10)
res.a_livestock <- get_model(path)

#' #### Posterior
#+ res.a_livestock, fig.height = 10
res.a_livestock %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_livestock %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_livestock %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_livestock %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_livestock %>% DIC() # 6151.664
res.a_livestock %>% testSSEF()
res.a_livestock %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_livestock %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_livestock %>% plot_autocorr()

#' #### Model
res.a_livestock$model

#' #### Results
res.a_livestock




#' ### res.a_cattle {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_cattle.R"),
                            thin = 10)
res.a_cattle <- get_model(path)

#' #### Posterior
#+ res.a_cattle, fig.height = 10
res.a_cattle %>% plot_density(data_ah)

#' #### Statistics
res.a_cattle %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_cattle %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_cattle %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_cattle %>% DIC() # 6151.664
res.a_cattle %>% testSSEF()
res.a_cattle %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_cattle %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_cattle %>% plot_autocorr()

#' #### Model
res.a_cattle$model

#' #### Results
res.a_cattle




#' ### res.a_pig {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_pig.R"), thin = 10)
res.a_pig <- get_model(path)

#' #### Posterior
#+ res.a_pig, fig.height = 10
res.a_pig %>% plot_density(data_ah)

#' #### Statistics
res.a_pig %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_pig %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_pig %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_pig %>% DIC() # 6151.664
res.a_pig %>% testSSEF()
res.a_pig %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_pig %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_pig %>% plot_autocorr()

#' #### Model
res.a_pig$model

#' #### Results
res.a_pig




#' ### res.a_chicken {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_chicken.R"),
                            thin = 10)
res.a_chicken <- get_model(path)

#' #### Posterior
#+ res.a_chicken, fig.height = 10
res.a_chicken %>% plot_density(data_ah)

#' #### Statistics
res.a_chicken %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_chicken %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_chicken %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_chicken %>% DIC() # 6151.664
res.a_chicken %>% testSSEF()
res.a_chicken %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_chicken %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_chicken %>% plot_autocorr()

#' #### Model
res.a_chicken$model

#' #### Results
res.a_chicken




#' ### res.a_livestock_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory,
                                               "a_livestock_subsets.R"),
                            thin = 10)
res.a_livestock_subsets <- get_model(path)

#' #### Posterior
#+ res.a_livestock_subsets, fig.height = 10
res.a_livestock_subsets %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_livestock_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_livestock_subsets %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_livestock_subsets %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_livestock_subsets %>% DIC() # 6151.664
res.a_livestock_subsets %>% testSSEF()
res.a_livestock_subsets %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_livestock_subsets %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_livestock_subsets %>% plot_autocorr()

#' #### Model
res.a_livestock_subsets$model

#' #### Results
res.a_livestock_subsets




#' ### res.a_companion {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_companion.R"),
                            thin = 10)
res.a_companion <- get_model(path)

#' #### Posterior
#+ a_companion, fig.height = 10
res.a_companion %>% plot_density(data_ah)

#' #### Statistics
res.a_companion %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_companion %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_companion %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_companion %>% DIC() # 6151.664
res.a_companion %>% testSSEF()
res.a_companion %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_companion %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_companion %>% plot_autocorr()

#' #### Model
res.a_companion$model

#' #### Results
res.a_companion




#' ### res.a_companion_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory,
                                               "a_companion_subsets.R"),
                            thin = 10)
res.a_companion_subsets <- get_model(path)

#' #### Posterior
#+ res.a_companion_subsets, fig.height = 10
res.a_companion_subsets %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_companion_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_companion_subsets %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_companion_subsets %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_companion_subsets %>% DIC() # 6151.664
res.a_companion_subsets %>% testSSEF()
res.a_companion_subsets %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_companion_subsets %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_companion_subsets %>% plot_autocorr()

#' #### Model
res.a_companion_subsets$model

#' #### Results
res.a_companion_subsets




#' ### res.a_wild {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_wild.R"),
                            thin = 10)
res.a_wild <- get_model(path)

#' #### Posterior
#+ a_wild, fig.height = 10
res.a_wild %>% plot_density(data_ah)

#' #### Statistics
res.a_wild %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_wild %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_wild %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_wild %>% DIC() # 6151.664
res.a_wild %>% testSSEF()
res.a_wild %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_wild %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_wild %>% plot_autocorr()

#' #### Model
res.a_wild$model

#' #### Results
res.a_wild




#' ### res.a_wild_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_wild_subsets.R"),
                            thin = 10)
res.a_wild_subsets <- get_model(path)

#' #### Posterior
#+ res.a_wild_subsets, fig.height = 10
res.a_wild_subsets %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_wild_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_wild_subsets %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_wild_subsets %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_wild_subsets %>% DIC() # 6151.664
res.a_wild_subsets %>% testSSEF()
res.a_wild_subsets %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_wild_subsets %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_wild_subsets %>% plot_autocorr()

#' #### Model
res.a_wild_subsets$model

#' #### Results
res.a_wild_subsets




#' ### res.a_types {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_types.R"),
                            thin = 10)
res.a_types <- get_model(path)

#' #### Posterior
#+ res.a_types, fig.height = 10
res.a_types %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_types %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_types %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_types %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_types %>% DIC() # 6151.664
res.a_types %>% testSSEF()
res.a_types %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_types %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_types %>% plot_autocorr()

#' #### Model
res.a_types$model

#' #### Results
res.a_types




#' ### res.a_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data_ah, file.path(directory, "a_subsets.R"),
                            thin = 10)
res.a_subsets <- get_model(path)

#' #### Posterior
#+ res.a_subsets, fig.height = 10
res.a_subsets %>% plot_density(data_ah)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_subsets %>% plot_antibiotics(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_subsets %>% plot_correlation(data_ah)

#' #### Diagnostics
res.a_subsets %>% DIC() # 6151.664
res.a_subsets %>% testSSEF()
res.a_subsets %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_subsets %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_subsets %>% plot_autocorr()

#' #### Model
res.a_subsets$model

#' #### Results
res.a_subsets




#' ### res.ac1 {.tabset}
#' response ~ antibiotic.class_{goodbad}
#' goodbad ~ clinical
#'

path <- run_SpARKjags_model(data, file.path(directory, "ac1.R"))
res.ac1 <- get_model(path)

#' #### Posterior
#+ res.ac1, fig.height = 10
res.ac1 %>% plot_density(data)

#' #### Statistics
res.ac1 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ac1 %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ac1 %>% plot_correlation(data)

#' #### Diagnostics
res.ac1 %>% DIC() # 6156.994
res.ac1 %>% testSSEF()
res.ac1 %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.ac1 %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.ac1 %>% plot_autocorr()

#' #### Model
res.ac1$model

#' #### Results
res.ac1




#' ### res.ac2 {.tabset}
#' response ~ antibiotic.class_{goodbad} + clinical
#'

path <- run_SpARKjags_model(data, file.path(directory, "ac2.R"))
res.ac2 <- get_model(path)

#' #### Posterior
#+ res.ac2, fig.height = 10
res.ac2 %>% plot_density(data)

#' #### Statistics
res.ac2 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ac2 %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ac2 %>% plot_correlation(data)

#' #### Diagnostics
res.ac2 %>% DIC() # 6079.197
res.ac2 %>% testSSEF()
res.ac2 %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.ac2 %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.ac2 %>% plot_autocorr()

#' #### Model
res.ac2$model

#' #### Results
res.ac2




#' ### res.a_c {.tabset}
#' response ~ antibiotic.class_{goodbad,clinical}
#'

path <- run_SpARKjags_model(data, file.path(directory, "a_c.R"))
res.a_c <- get_model(path)

#' #### Posterior
#+ res.a_c, fig.height = 10
res.a_c %>% plot_density(data)

#' #### Statistics
res.a_c %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.a_c %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.a_c %>% plot_correlation(data)

#' #### Diagnostics
res.a_c %>% DIC() # 6045.375
res.a_c %>% testSSEF()
res.a_c %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.a_c %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.a_c %>% plot_autocorr()

#' #### Model
res.a_c$model

#' #### Results
res.a_c


#' ### Summary of DIC results {-}
#'
DICtable(c("res.ac1", "res.ac2", "res.a_c"))




#' ### res.asm {.tabset}
#' response ~ antibiotic.class_{goodbad} + sample.month
#'

path <- run_SpARKjags_model(data, file.path(directory, "asm.R"), thin = 10)
res.asm <- get_model(path)

#' #### Posterior
#+ res.asm, fig.height = 10
res.asm %>% plot_density(data)

#' #### Statistics
res.asm %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.asm %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.asm %>% plot_correlation(data)

#' #### Diagnostics
res.asm %>% DIC() # 6065.661
res.asm %>% testSSEF()
res.asm %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.asm %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.asm %>% plot_autocorr()

#' #### Model
res.asm$model

#' #### Results
res.asm




#' ### res.ass {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season
#'

path <- run_SpARKjags_model(data, file.path(directory, "ass.R"), thin = 30)
res.ass <- get_model(path)

#' #### Posterior
#+ res.ass, fig.height = 10
res.ass %>% plot_density(data)

#' #### Statistics
res.ass %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ass %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ass %>% plot_correlation(data)

#' #### Diagnostics
res.ass %>% DIC() # 6049.671
res.ass %>% testSSEF()
res.ass %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.ass %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.ass %>% plot_autocorr()

#' #### Model
res.ass$model

#' #### Results
res.ass


#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))




#'
#' ## 1b. Demographic
#' Find the best model with gender and age
#'
#' ### res.assg {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + gender
#'

path <- run_SpARKjags_model(data, file.path(directory, "assg.R"), thin = 30)
res.assg <- get_model(path)

#' #### Posterior
#+ res.assg, fig.height = 10
res.assg %>% plot_density(data)

#' #### Statistics
res.assg %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assg %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assg %>% plot_correlation(data)

#' #### Diagnostics
res.assg %>% DIC() # 6035.019
res.assg %>% testSSEF()
res.assg %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assg %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assg %>% plot_autocorr()

#' #### Model
res.assg$model

#' #### Results
res.assg




#' ### res.assag {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup
#'

path <- run_SpARKjags_model(data, file.path(directory, "assag.R"), thin = 30)
res.assag <- get_model(path)

#' #### Posterior
#+ res.assag, fig.height = 10
res.assag %>% plot_density(data)

#' #### Statistics
res.assag %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assag %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assag %>% plot_correlation(data)

#' #### Diagnostics
res.assag %>% DIC() # 6030.135
res.assag %>% testSSEF()
res.assag %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assag %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assag %>% plot_autocorr()

#' #### Model
res.assag$model

#' #### Results
res.assag




#' ### res.assag2 {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup2
#'

path <- run_SpARKjags_model(data, file.path(directory, "assag2.R"), thin = 30)
res.assag2 <- get_model(path)

#' #### Posterior
#+ res.assag2, fig.height = 10
res.assag2 %>% plot_density(data)

#' #### Statistics
res.assag2 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assag2 %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assag2 %>% plot_correlation(data)

#' #### Diagnostics
res.assag2 %>% DIC() # 6047.175
res.assag2 %>% testSSEF()
res.assag2 %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assag2 %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assag2 %>% plot_autocorr()

#' #### Model
res.assag2$model

#' #### Results
res.assag2




#' ### res.assage {.tabset}
#' response ~ antibiotic_class_{goodbad} + sampling_month + age
#'

path <- run_SpARKjags_model(data, file.path(directory, "assage.R"), thin = 20)
res.assage <- get_model(path)

#' #### Posterior
#+ res.assage, fig.height = 10
res.assage %>% plot_density(data)

#' #### Statistics
res.assage %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assage %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assage %>% plot_correlation(data)

#' #### Diagnostics
res.assage %>% DIC() # 6032.315
res.assage %>% testSSEF()
res.assage %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assage %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assage %>% plot_autocorr()

#' #### Model
res.assage$model

#' #### Results
res.assage




#' ### res.assagesq {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + age^2 + age
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagesq.R"), thin = 20)
res.assagesq <- get_model(path)

#' #### Posterior
#+ res.assagesq, fig.height = 10
res.assagesq %>% plot_density(data)

#' #### Statistics
res.assagesq %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagesq %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagesq %>% plot_correlation(data)

#' #### Diagnostics
res.assagesq %>% DIC() # 6095.407
res.assagesq %>% testSSEF()
res.assagesq %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagesq %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagesq %>% plot_autocorr()

#' #### Model
res.assagesq$model

#' #### Results
res.assagesq





#' ### res.assagg {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + gender
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagg.R"), thin = 20)
res.assagg <- get_model(path)

#' #### Posterior
#+ res.assagg, fig.height = 10
res.assagg %>% plot_density(data)

#' #### Statistics
res.assagg %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagg %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagg %>% plot_correlation(data)

#' #### Diagnostics
res.assagg %>% DIC() # 6049.079
res.assagg %>% testSSEF()
res.assagg %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagg %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagg %>% plot_autocorr()

#' #### Model
res.assagg$model

#' #### Results
res.assagg


#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))
DICtable(c("res.assg", "res.assag", "res.assag2", "res.assage",
           "res.assagesq", "res.assagg"))




#'
#' ## 1c. Hospital structure
#' Find the best model with hospital, wardtype, and ward
#'
#' ### res.assagh {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + hospital
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagh.R"), thin = 20)
res.assagh <- get_model(path)

#' #### Posterior
#+ res.assagh, fig.height = 10
res.assagh %>% plot_density(data)

#' #### Statistics
res.assagh %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagh %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagh %>% plot_correlation(data)

#' #### Diagnostics
res.assagh %>% DIC() # 5998.832
res.assagh %>% testSSEF()
res.assagh %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagh %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagh %>% plot_autocorr()

#' #### Model
res.assagh$model

#' #### Results
res.assagh




#' ### res.assagwt {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + wardtype
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagwt.R"), thin = 20)
res.assagwt <- get_model(path)

#' #### Posterior
#+ res.assagwt, fig.height = 10
res.assagwt %>% plot_density(data)

#' #### Statistics
res.assagwt %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwt %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwt %>% plot_correlation(data)

#' #### Diagnostics
res.assagwt %>% DIC() # 5953.609
res.assagwt %>% testSSEF()
res.assagwt %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagwt %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagwt %>% plot_autocorr()

#' #### Model
res.assagwt$model

#' #### Results
res.assagwt




#' ### res.assagw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagw.R"), thin = 20)
res.assagw <- get_model(path)

#' #### Posterior
#+ res.assagw, fig.height = 10
res.assagw %>% plot_density(data)

#' #### Statistics
res.assagw %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagw %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagw %>% plot_correlation(data)

#' #### Diagnostics
res.assagw %>% DIC() # 5883.045
res.assagw %>% testSSEF()
res.assagw %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagw %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagw %>% plot_autocorr()

#' #### Model
res.assagw$model

#' #### Results
res.assagw




#' ### res.assagwt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup +
#' wardtype_{ward}
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagwt_w.R"),
                            thin = 20)
res.assagwt_w <- get_model(path)

#' #### Posterior
#+ res.assagwt_w, fig.height = 10
res.assagwt_w %>% plot_density(data)

#' #### Statistics
res.assagwt_w %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwt_w %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwt_w %>% plot_correlation(data)

#' #### Diagnostics
res.assagwt_w %>% DIC() # 5882.934
res.assagwt_w %>% testSSEF()
res.assagwt_w %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagwt_w %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagwt_w %>% plot_autocorr()

#' #### Model
res.assagwt_w$model

#' #### Results
res.assagwt_w



#' ### res.assaghwtw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + hospital +
#' wardtype + ward
#'

path <- run_SpARKjags_model(data, file.path(directory, "assaghwtw.R"),
                            thin = 20)
res.assaghwtw <- get_model(path)

#' #### Posterior
#+ res.assaghwtw, fig.height = 10
res.assaghwtw %>% plot_density(data)

#' #### Statistics
res.assaghwtw %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assaghwtw %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assaghwtw %>% plot_correlation(data)

#' #### Diagnostics
res.assaghwtw %>% DIC() # 5885.554
res.assaghwtw %>% testSSEF()
res.assaghwtw %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assaghwtw %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assaghwtw %>% plot_autocorr()

#' #### Model
res.assaghwtw$model

#' #### Results
res.assaghwtw




#' ### res.assagh_wt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup +
#' hospital_{wardtype_{ward}}
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagh_wt_w.R"),
                            thin = 20)
res.assagh_wt_w <- get_model(path)

#' #### Posterior
#+ res.assagh_wt_w, fig.height = 10
res.assagh_wt_w %>% plot_density(data)

#' #### Statistics
res.assagh_wt_w %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagh_wt_w %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagh_wt_w %>% plot_correlation(data)

#' #### Diagnostics
res.assagh_wt_w %>% DIC() # 5882.687
res.assagh_wt_w %>% testSSEF()
res.assagh_wt_w %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagh_wt_w %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagh_wt_w %>% plot_autocorr()

#' #### Model
res.assagh_wt_w$model

#' #### Results
res.assagh_wt_w


#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))
DICtable(c("res.assg", "res.assag", "res.assag2", "res.assage",
           "res.assagesq", "res.assagg"))
DICtable(c("res.assagh", "res.assagwt", "res.assagw", "res.assagwt_w",
           "res.assaghwtw", "res.assagh_wt_w"))

#' Choose res.assagw
#'




#' ### res.assagwc {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagwc.R"), thin = 20)
res.assagwc <- get_model(path)

#' #### Posterior
#+ res.assagwc, fig.height = 10
res.assagwc %>% plot_density(data)

#' #### Statistics
res.assagwc %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwc %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwc %>% plot_correlation(data)

#' #### Diagnostics
res.assagwc %>% DIC() # 5886.387
res.assagwc %>% testSSEF()
res.assagwc %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagwc %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagwc %>% plot_autocorr()

#' #### Model
res.assagwc$model

#' #### Results
res.assagwc




#' ### res.assagwcst {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical + sample_type
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagwcst.R"), thin = 30)
res.assagwcst <- get_model(path)

#' #### Posterior
#+ res.assagwcst, fig.height = 10
res.assagwcst %>% plot_density(data)

#' #### Statistics
res.assagwcst %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwcst %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwcst %>% plot_correlation(data)

#' #### Diagnostics
res.assagwcst %>% DIC() # 5835.189
res.assagwcst %>% testSSEF()
res.assagwcst %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagwcst %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagwcst %>% plot_autocorr()

#' #### Model
res.assagwcst$model

#' #### Results
res.assagwcst




#' ### res.assagwc_st {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical_{sampletype}
#'

path <- run_SpARKjags_model(data, file.path(directory, "assagwc_st.R"),
                            thin = 30)
res.assagwc_st <- get_model(path)

#' #### Posterior
#+ res.assagwc_st, fig.height = 10
res.assagwc_st %>% plot_density(data)

#' #### Statistics
res.assagwc_st %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwc_st %>% plot_antibiotics(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwc_st %>% plot_correlation(data)

#' #### Diagnostics
res.assagwc_st %>% DIC() # 5846.981
res.assagwc_st %>% testSSEF()
res.assagwc_st %>% testPSRF()

#' #### Trace plot
#+ fig.height = 6
res.assagwc_st %>% plot_caterpillar()

#' #### Autocorrelation
#+ fig.height = 6
res.assagwc_st %>% plot_autocorr()

#' #### Model
res.assagwc_st$model

#' #### Results
res.assagwc_st



#' #' ### res.a_cssagwst {.tabset}
#' #' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup +
#' #' ward + sampletype
#' #'
#'
#' path <- run_SpARKjags_model(data, file.path(directory, "a_cssagwst.R"),
#'                             thin = 10)
#' res.a_cssagwst <- get_model(path)
#'
#' #' #### Posterior
#' #+ res.a_cssagwst, fig.height = 10
#' res.a_cssagwst %>% plot_density(data)
#'
#' #' #### Statistics
#' res.a_cssagwst %>% summarise_samples(data)
#'
#' #' #### AMR Summary
#' #+ fig.height = 10
#' res.a_cssagwst %>% plot_antibiotics(data)
#'
#' #' #### AMR Correlation
#' #+ fig.height = 10
#' res.a_cssagwst %>% plot_correlation(data)
#'
#' #' #### Diagnostics
#' res.a_cssagwst %>% DIC() # 5773.768
#' res.a_cssagwst %>% testSSEF()
#' res.a_cssagwst %>% testPSRF()
#'
#' #' #### Trace plot
#' #+ fig.height = 6
#' res.a_cssagwst %>% plot_caterpillar()
#'
#' #' #### Autocorrelation
#' #+ fig.height = 6
#' res.a_cssagwst %>% plot_autocorr()
#'
#' #' #### Model
#' res.a_cssagwst$model
#'
#' #' #### Results
#' res.a_cssagwst


#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))
DICtable(c("res.assg", "res.assag", "res.assag2", "res.assage",
           "res.assagesq", "res.assagg"))
DICtable(c("res.assagh", "res.assagwt", "res.assagw", "res.assagwt_w",
           "res.assaghwtw", "res.assagh_wt_w"))
DICtable(c("res.assagwc", "res.assagwcst", "res.assagwc_st"))
# "res.a_cssagwst"

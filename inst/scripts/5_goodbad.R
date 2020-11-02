#' ---
#' title: Goodbad models
#' date: Feb 2020
#' output:
#'   html_document:
#'     theme: paper
#'     format: html_clean
#'     code_folding: hide
#'     highlight: pygments
#'     toc: TRUE
#'     toc_float: TRUE
#' ---
#' Last updated: `r Sys.time()`

#+ setup, include=F

library(SpARKjags)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE)

# Find models here (installed with the package)
directory <- "goodbad_models"
# Save results here (make sure Desktop/git/SpARK/SpARKjags is your working
# directory)
res_dir <- file.path(getwd(), "results", directory)

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
#' * Posterior: The tails of the violins are trimmed to the range of the data.
#'   All violins are scaled to have the same maximum width. The grey line
#'   represents the 1st and 3rd quantile distance around the median
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

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a_naive.R"),
                            save_to = res_dir)
res.a_naive <- get_model(path)
model_name <- "res.a_naive"

#' #### Posterior
#+ res.a_naive, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_naive$model

#' #### Results
res.a_naive




#' ### res.a {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a.R"),
                            save_to = res_dir,
                            thin = 10)
res.a <- get_model(path)
model_name <- "res.a"

#' #### Posterior
#+ res.a, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_livestock.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_livestock <- get_model(path)
model_name <- "res.a_livestock"

#' #### Posterior
#+ res.a_livestock, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_livestock$model

#' #### Results
res.a_livestock




#' ### res.a_cattle {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_cattle.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_cattle <- get_model(path)
model_name <- "res.a_cattle"

#' #### Posterior
#+ res.a_cattle, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_cattle$model

#' #### Results
res.a_cattle




#' ### res.a_pig {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_pig.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_pig <- get_model(path)
model_name <- "res.a_pig"

#' #### Posterior
#+ res.a_pig, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_pig$model

#' #### Results
res.a_pig




#' ### res.a_chicken {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_chicken.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_chicken <- get_model(path)
model_name <- "res.a_chicken"

#' #### Posterior
#+ res.a_chicken, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_chicken$model

#' #### Results
res.a_chicken




#' ### res.a_livestock_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory,
                                                        "a_livestock_subsets.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_livestock_subsets <- get_model(path)
model_name <- "res.a_livestock_subsets"

#' #### Posterior
#+ res.a_livestock_subsets, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_livestock_subsets$model

#' #### Results
res.a_livestock_subsets




#' ### res.a_companion {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_companion.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_companion <- get_model(path)
model_name <- "res.a_companion"

#' #### Posterior
#+ a_companion, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_companion$model

#' #### Results
res.a_companion




#' ### res.a_companion_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory,
                                                        "a_companion_subsets.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_companion_subsets <- get_model(path)
model_name <- "res.a_companion_subsets"

#' #### Posterior
#+ res.a_companion_subsets, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_companion_subsets$model

#' #### Results
res.a_companion_subsets




#' ### res.a_wild {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_wild.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_wild <- get_model(path)
model_name <- "res.a_wild"

#' #### Posterior
#+ a_wild, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_wild$model

#' #### Results
res.a_wild




#' ### res.a_wild_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_wild_subsets.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_wild_subsets <- get_model(path)
model_name <- "res.a_wild_subsets"

#' #### Posterior
#+ res.a_wild_subsets, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_wild_subsets$model

#' #### Results
res.a_wild_subsets




#' ### res.a_types {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_types.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_types <- get_model(path)
model_name <- "res.a_types"

#' #### Posterior
#+ res.a_types, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_types$model

#' #### Results
res.a_types




#' ### res.a_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory, "a_subsets.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_subsets <- get_model(path)
model_name <- "res.a_subsets"

#' #### Posterior
#+ res.a_subsets, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data_ah)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data_ah)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.a_subsets$model

#' #### Results
res.a_subsets




#' ### res.ac1 {.tabset}
#' response ~ antibiotic.class_{goodbad}
#' goodbad ~ clinical
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ac1.R"),
                            save_to = res_dir)
res.ac1 <- get_model(path)
model_name <- "res.ac1"

#' #### Posterior
#+ res.ac1, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.ac1$model

#' #### Results
res.ac1




#' ### res.ac2 {.tabset}
#' response ~ antibiotic.class_{goodbad} + clinical
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ac2.R"),
                            save_to = res_dir)
res.ac2 <- get_model(path)
model_name <- "res.ac2"

#' #### Posterior
#+ res.ac2, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.ac2$model

#' #### Results
res.ac2




#' ### res.a_c {.tabset}
#' response ~ antibiotic.class_{goodbad,clinical}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a_c.R"),
                            save_to = res_dir)
res.a_c <- get_model(path)
model_name <- "res.a_c"

#' #### Posterior
#+ res.a_c, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "asm.R"),
                            save_to = res_dir,
                            thin = 10)
res.asm <- get_model(path)
model_name <- "res.asm"

#' #### Posterior
#+ res.asm, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.asm$model

#' #### Results
res.asm




#' ### res.ass {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ass.R"),
                            save_to = res_dir,
                            thin = 30)
res.ass <- get_model(path)
model_name <- "res.ass"

#' #### Posterior
#+ res.ass, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assg.R"),
                            save_to = res_dir,
                            thin = 30)
res.assg <- get_model(path)
model_name <- "res.assg"

#' #### Posterior
#+ res.assg, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assg$model

#' #### Results
res.assg




#' ### res.assag {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assag.R"),
                            save_to = res_dir,
                            thin = 30)
res.assag <- get_model(path)
model_name <- "res.assag"

#' #### Posterior
#+ res.assag, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assag$model

#' #### Results
res.assag




#' ### res.assag2 {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup2
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assag2.R"),
                            save_to = res_dir,
                            thin = 30)
res.assag2 <- get_model(path)
model_name <- "res.assag2"

#' #### Posterior
#+ res.assag2, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assag2$model

#' #### Results
res.assag2




#' ### res.assage {.tabset}
#' response ~ antibiotic_class_{goodbad} + sampling_month + age
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assage.R"),
                            save_to = res_dir,
                            thin = 20)
res.assage <- get_model(path)
model_name <- "res.assage"

#' #### Posterior
#+ res.assage, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assage$model

#' #### Results
res.assage




#' ### res.assagesq {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + age^2 + age
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagesq.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagesq <- get_model(path)
model_name <- "res.assagesq"

#' #### Posterior
#+ res.assagesq, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagesq$model

#' #### Results
res.assagesq





#' ### res.assagg {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + gender
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagg.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagg <- get_model(path)
model_name <- "res.assagg"

#' #### Posterior
#+ res.assagg, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagh.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagh <- get_model(path)
model_name <- "res.assagh"

#' #### Posterior
#+ res.assagh, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))
#' #### Model
res.assagh$model

#' #### Results
res.assagh




#' ### res.assagwt {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + wardtype
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagwt.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagwt <- get_model(path)
model_name <- "res.assagwt"

#' #### Posterior
#+ res.assagwt, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagwt$model

#' #### Results
res.assagwt




#' ### res.assagw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagw.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagw <- get_model(path)
model_name <- "res.assagw"

#' #### Posterior
#+ res.assagw, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagw$model

#' #### Results
res.assagw




#' ### res.assagwt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup +
#' wardtype_{ward}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagwt_w.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagwt_w <- get_model(path)
model_name <- "res.assagwt_w"

#' #### Posterior
#+ res.assagwt_w, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagwt_w$model

#' #### Results
res.assagwt_w



#' ### res.assaghwtw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + hospital +
#' wardtype + ward
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assaghwtw.R"),
                            save_to = res_dir,
                            thin = 20)
res.assaghwtw <- get_model(path)
model_name <- "res.assaghwtw"

#' #### Posterior
#+ res.assaghwtw, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assaghwtw$model

#' #### Results
res.assaghwtw




#' ### res.assagh_wt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup +
#' hospital_{wardtype_{ward}}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagh_wt_w.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagh_wt_w <- get_model(path)
model_name <- "res.assagh_wt_w"

#' #### Posterior
#+ res.assagh_wt_w, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagwc.R"),
                            save_to = res_dir,
                            thin = 20)
res.assagwc <- get_model(path)
model_name <- "res.assagwc"

#' #### Posterior
#+ res.assagwc, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagwc$model

#' #### Results
res.assagwc




#' ### res.assagwcst {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical + sample_type
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagwcst.R"),
                            save_to = res_dir,
                            thin = 30)
res.assagwcst <- get_model(path)
model_name <- "res.assagwcst"

#' #### Posterior
#+ res.assagwcst, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
res.assagwcst$model

#' #### Results
res.assagwcst




#' ### res.assagwc_st {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical_{sampletype}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "assagwc_st.R"),
                            save_to = res_dir,
                            thin = 30)
res.assagwc_st <- get_model(path)
model_name <- "res.assagwc_st"

#' #### Posterior
#+ res.assagwc_st, fig.height = 10
plot_density(model = get(model_name),
             data = data,
             save_to = res_dir,
             model_name = model_name)

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 10
plot_antibiotics(model = get(model_name),
                 data = data)

#' #### AMR Correlation
#+ fig.height = 10
plot_correlation(model = get(model_name),
                 data = data)

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = res_dir)

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

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

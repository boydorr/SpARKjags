#' ---
#' title: Seemingly unrelated regressions
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
library(SpARK)
library(SpARKcarbapenem)
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
#' * Posterior: The grey bar represents the 1st and 3rd quantile distance around the median
#' * AMR Correlation: calculated Pearson correlation using pairwise.complete.obs
#' * Class interpretation: R when sample is resistant to at least one antibiotic, S when sample is susceptible to at least one antibiotic, otherwise NA
#' * Samples are determined as being in the bad group when bad.p (bad.gp, bad,v, and bad.o) is greater than 0.5
#'
#' ## 1a. Experimental structure
#'
#' ### res.a_naive {.tabset}
#' response ~ antibiotic_class (naive model)
#'


# run_model(data, "a_naive", directory)
res.a_naive <- get_model("a_naive", directory)

#' #### Posterior
#+ res.a_naive, fig.height = 6
var.regex <- "(a.prob)|(intercept)|(sd)"
params <- list(`probability of resistance` = "prob", "intercept", "sd")
tmp <- data$lookup$antibiotic_class %>%
  dplyr::mutate(index = paste0("a.prob[", 1:13, "]"))
labels <- list(tmp, NA, NA)
res.a_naive %>% densityplot(data, var.regex, params, labels)

#' #### Diagnostics
res.a_naive %>% DIC() # 6166.735
res.a_naive %>% testSSEF()
res.a_naive %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_naive %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_naive %>% plotAC(var.regex)

#' #### Model
res.a_naive$model

#' #### Results
res.a_naive




#' ### res.a {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data, "a", directory, thin = 10)
res.a <- get_model("a", directory)

#' #### Posterior
#+ res.a, fig.height = 10
var.regex <- get_vars(res.a)
params <- get_params()
labels <- get_labels(data)
res.a %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.a %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.a %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.a %>% plotcorrelation(data)

#' #### Diagnostics
res.a %>% DIC() # 6151.664
res.a %>% testSSEF()
res.a %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a %>% plotAC(var.regex)

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
  ggplot2::geom_bar(ggplot2::aes(x = ASSOCIATED_SPECIES, y = count, fill = TYPE),
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

# run_model(data_ah, "a_livestock", directory, thin = 10)
res.a_livestock <- get_model("a_livestock", directory)

#' #### Posterior
#+ res.a_livestock, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_livestock %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_livestock %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_livestock %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_livestock %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_livestock %>% DIC() # 6151.664
res.a_livestock %>% testSSEF()
res.a_livestock %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_livestock %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_livestock %>% plotAC(var.regex)

#' #### Model
res.a_livestock$model

#' #### Results
res.a_livestock




#' ### res.a_cattle {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_cattle", directory, thin = 10)
res.a_cattle <- get_model("a_cattle", directory)

#' #### Posterior
#+ res.a_cattle, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_cattle %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
res.a_cattle %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_cattle %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_cattle %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_cattle %>% DIC() # 6151.664
res.a_cattle %>% testSSEF()
res.a_cattle %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_cattle %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_cattle %>% plotAC(var.regex)

#' #### Model
res.a_cattle$model

#' #### Results
res.a_cattle




#' ### res.a_pig {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_pig", directory, thin = 10)
res.a_pig <- get_model("a_pig", directory)

#' #### Posterior
#+ res.a_pig, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_pig %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
res.a_pig %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_pig %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_pig %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_pig %>% DIC() # 6151.664
res.a_pig %>% testSSEF()
res.a_pig %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_pig %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_pig %>% plotAC(var.regex)

#' #### Model
res.a_pig$model

#' #### Results
res.a_pig




#' ### res.a_chicken {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_chicken", directory, thin = 10)
res.a_chicken <- get_model("a_chicken", directory)

#' #### Posterior
#+ res.a_chicken, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_chicken %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
res.a_chicken %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_chicken %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_chicken %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_chicken %>% DIC() # 6151.664
res.a_chicken %>% testSSEF()
res.a_chicken %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_chicken %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_chicken %>% plotAC(var.regex)

#' #### Model
res.a_chicken$model

#' #### Results
res.a_chicken




#' ### res.a_livestock_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_livestock_subsets", directory, thin = 10)
res.a_livestock_subsets <- get_model("a_livestock_subsets", directory)

#' #### Posterior
#+ res.a_livestock_subsets, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_livestock_subsets %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_livestock_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_livestock_subsets %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_livestock_subsets %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_livestock_subsets %>% DIC() # 6151.664
res.a_livestock_subsets %>% testSSEF()
res.a_livestock_subsets %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_livestock_subsets %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_livestock_subsets %>% plotAC(var.regex)

#' #### Model
res.a_livestock_subsets$model

#' #### Results
res.a_livestock_subsets




#' ### res.a_companion {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_companion", directory, thin = 10)
res.a_companion <- get_model("a_companion", directory)

#' #### Posterior
#+ a_companion, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_companion %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
res.a_companion %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_companion %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_companion %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_companion %>% DIC() # 6151.664
res.a_companion %>% testSSEF()
res.a_companion %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_companion %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_companion %>% plotAC(var.regex)

#' #### Model
res.a_companion$model

#' #### Results
res.a_companion




#' ### res.a_companion_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_companion_subsets", directory, thin = 10)
res.a_companion_subsets <- get_model("a_companion_subsets", directory)

#' #### Posterior
#+ res.a_companion_subsets, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_companion_subsets %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_companion_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_companion_subsets %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_companion_subsets %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_companion_subsets %>% DIC() # 6151.664
res.a_companion_subsets %>% testSSEF()
res.a_companion_subsets %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_companion_subsets %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_companion_subsets %>% plotAC(var.regex)

#' #### Model
res.a_companion_subsets$model

#' #### Results
res.a_companion_subsets




#' ### res.a_wild {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_wild", directory, thin = 10)
res.a_wild <- get_model("a_wild", directory)

#' #### Posterior
#+ a_wild, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_wild %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
res.a_wild %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_wild %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_wild %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_wild %>% DIC() # 6151.664
res.a_wild %>% testSSEF()
res.a_wild %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_wild %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_wild %>% plotAC(var.regex)

#' #### Model
res.a_wild$model

#' #### Results
res.a_wild




#' ### res.a_wild_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_wild_subsets", directory, thin = 10)
res.a_wild_subsets <- get_model("a_wild_subsets", directory)

#' #### Posterior
#+ res.a_wild_subsets, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_wild_subsets %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_wild_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_wild_subsets %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_wild_subsets %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_wild_subsets %>% DIC() # 6151.664
res.a_wild_subsets %>% testSSEF()
res.a_wild_subsets %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_wild_subsets %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_wild_subsets %>% plotAC(var.regex)

#' #### Model
res.a_wild_subsets$model

#' #### Results
res.a_wild_subsets




#' ### res.a_types {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_types", directory, thin = 10)
res.a_types <- get_model("a_types", directory)

#' #### Posterior
#+ res.a_types, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_types %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_types %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_types %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_types %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_types %>% DIC() # 6151.664
res.a_types %>% testSSEF()
res.a_types %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_types %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_types %>% plotAC(var.regex)

#' #### Model
res.a_types$model

#' #### Results
res.a_types




#' ### res.a_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

# run_model(data_ah, "a_subsets", directory, thin = 10)
res.a_subsets <- get_model("a_subsets", directory)

#' #### Posterior
#+ res.a_subsets, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data_ah)
res.a_subsets %>% densityplot(data_ah, var.regex, params, labels)

#' #### Statistics
# Posterior probability of each sample being in the bad group
res.a_subsets %>% summarise_samples(data_ah)

#' #### AMR Summary
#+ fig.height = 10
res.a_subsets %>% antibioticsplot(data_ah)

#' #### AMR Correlation
#+ fig.height = 10
res.a_subsets %>% plotcorrelation(data_ah)

#' #### Diagnostics
res.a_subsets %>% DIC() # 6151.664
res.a_subsets %>% testSSEF()
res.a_subsets %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_subsets %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_subsets %>% plotAC(var.regex)

#' #### Model
res.a_subsets$model

#' #### Results
res.a_subsets




#' ### res.ac1 {.tabset}
#' response ~ antibiotic.class_{goodbad}
#' goodbad ~ clinical
#'

# run_model(data, "ac1", directory)
res.ac1 <- get_model("ac1", directory)

#' #### Posterior
#+ res.ac1, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.ac1 %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.ac1 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ac1 %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ac1 %>% plotcorrelation(data)

#' #### Diagnostics
res.ac1 %>% DIC() # 6156.994
res.ac1 %>% testSSEF()
res.ac1 %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.ac1 %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.ac1 %>% plotAC(var.regex)

#' #### Model
res.ac1$model

#' #### Results
res.ac1




#' ### res.ac2 {.tabset}
#' response ~ antibiotic.class_{goodbad} + clinical
#'

# run_model(data, "ac2", directory)
res.ac2 <- get_model("ac2", directory)

#' #### Posterior
#+ res.ac2, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.ac2 %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.ac2 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ac2 %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ac2 %>% plotcorrelation(data)

#' #### Diagnostics
res.ac2 %>% DIC() # 6079.197
res.ac2 %>% testSSEF()
res.ac2 %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.ac2 %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.ac2 %>% plotAC(var.regex)

#' #### Model
res.ac2$model

#' #### Results
res.ac2




#' ### res.a_c {.tabset}
#' response ~ antibiotic.class_{goodbad,clinical}
#'

# run_model(data, "a_c", directory)
res.a_c <- get_model("a_c", directory)

#' #### Posterior
#+ res.a_c, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.a_c %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.a_c %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.a_c %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.a_c %>% plotcorrelation(data)

#' #### Diagnostics
res.a_c %>% DIC() # 6045.375
res.a_c %>% testSSEF()
res.a_c %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_c %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_c %>% plotAC(var.regex)

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

# run_model(data, "asm", directory, thin = 10)
res.asm <- get_model("asm", directory)

#' #### Posterior
#+ res.asm, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.asm %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.asm %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.asm %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.asm %>% plotcorrelation(data)

#' #### Diagnostics
res.asm %>% DIC() # 6065.661
res.asm %>% testSSEF()
res.asm %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.asm %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.asm %>% plotAC(var.regex)

#' #### Model
res.asm$model

#' #### Results
res.asm




#' ### res.ass {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season
#'

# run_model(data, "ass", directory, thin = 30)
res.ass <- get_model("ass", directory)

#' #### Posterior
#+ res.ass, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.ass %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.ass %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.ass %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.ass %>% plotcorrelation(data)

#' #### Diagnostics
res.ass %>% DIC() # 6049.671
res.ass %>% testSSEF()
res.ass %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.ass %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.ass %>% plotAC(var.regex)

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

# run_model(data, "assg", directory, thin = 30)
res.assg <- get_model("assg", directory)

#' #### Posterior
#+ res.assg, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assg %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assg %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assg %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assg %>% plotcorrelation(data)

#' #### Diagnostics
res.assg %>% DIC() # 6035.019
res.assg %>% testSSEF()
res.assg %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assg %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assg %>% plotAC(var.regex)

#' #### Model
res.assg$model

#' #### Results
res.assg




#' ### res.assag {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup
#'

# run_model(data, "assag", directory, thin = 30)
res.assag <- get_model("assag", directory)

#' #### Posterior
#+ res.assag, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assag %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assag %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assag %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assag %>% plotcorrelation(data)

#' #### Diagnostics
res.assag %>% DIC() # 6030.135
res.assag %>% testSSEF()
res.assag %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assag %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assag %>% plotAC(var.regex)

#' #### Model
res.assag$model

#' #### Results
res.assag




#' ### res.assag2 {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup2
#'

# run_model(data, "assag2", directory, thin = 30)
res.assag2 <- get_model("assag2", directory)

#' #### Posterior
#+ res.assag2, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assag2 %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assag2 %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assag2 %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assag2 %>% plotcorrelation(data)

#' #### Diagnostics
res.assag2 %>% DIC() # 6047.175
res.assag2 %>% testSSEF()
res.assag2 %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assag2 %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assag2 %>% plotAC(var.regex)

#' #### Model
res.assag2$model

#' #### Results
res.assag2




#' ### res.assage {.tabset}
#' response ~ antibiotic_class_{goodbad} + sampling_month + age
#'

# run_model(data, "assage", directory, thin = 20)
res.assage <- get_model("assage", directory)

#' #### Posterior
#+ res.assage, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assage %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assage %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assage %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assage %>% plotcorrelation(data)

#' #### Diagnostics
res.assage %>% DIC() # 6032.315
res.assage %>% testSSEF()
res.assage %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assage %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assage %>% plotAC(var.regex)

#' #### Model
res.assage$model

#' #### Results
res.assage




#' ### res.assagesq {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + age^2 + age
#'

# run_model(data, "assagesq", directory, thin = 20)
res.assagesq <- get_model("assagesq", directory)

#' #### Posterior
#+ res.assagesq, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagesq %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagesq %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagesq %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagesq %>% plotcorrelation(data)

#' #### Diagnostics
res.assagesq %>% DIC() # 6095.407
res.assagesq %>% testSSEF()
res.assagesq %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagesq %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagesq %>% plotAC(var.regex)

#' #### Model
res.assagesq$model

#' #### Results
res.assagesq





#' ### res.assagg {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + gender
#'

# run_model(data, "assagg", directory, thin = 20)
res.assagg <- get_model("assagg", directory)

#' #### Posterior
#+ res.assagg, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagg %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagg %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagg %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagg %>% plotcorrelation(data)

#' #### Diagnostics
res.assagg %>% DIC() # 6049.079
res.assagg %>% testSSEF()
res.assagg %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagg %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagg %>% plotAC(var.regex)

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

# run_model(data, "assagh", directory, thin = 20)
res.assagh <- get_model("assagh", directory)

#' #### Posterior
#+ res.assagh, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagh %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagh %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagh %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagh %>% plotcorrelation(data)

#' #### Diagnostics
res.assagh %>% DIC() # 5998.832
res.assagh %>% testSSEF()
res.assagh %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagh %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagh %>% plotAC(var.regex)

#' #### Model
res.assagh$model

#' #### Results
res.assagh




#' ### res.assagwt {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + wardtype
#'

# run_model(data, "assagwt", directory, thin = 20)
res.assagwt <- get_model("assagwt", directory)

#' #### Posterior
#+ res.assagwt, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagwt %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagwt %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwt %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwt %>% plotcorrelation(data)

#' #### Diagnostics
res.assagwt %>% DIC() # 5953.609
res.assagwt %>% testSSEF()
res.assagwt %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagwt %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagwt %>% plotAC(var.regex)

#' #### Model
res.assagwt$model

#' #### Results
res.assagwt




#' ### res.assagw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward
#'

# run_model(data, "assagw", directory, thin = 20)
res.assagw <- get_model("assagw", directory)

#' #### Posterior
#+ res.assagw, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagw %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagw %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagw %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagw %>% plotcorrelation(data)

#' #### Diagnostics
res.assagw %>% DIC() # 5883.045
res.assagw %>% testSSEF()
res.assagw %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagw %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagw %>% plotAC(var.regex)

#' #### Model
res.assagw$model

#' #### Results
res.assagw




#' ### res.assagwt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + wardtype_{ward}
#'

# run_model(data, "assagwt_w", directory, thin = 20)
res.assagwt_w <- get_model("assagwt_w", directory)

#' #### Posterior
#+ res.assagwt_w, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagwt_w %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagwt_w %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwt_w %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwt_w %>% plotcorrelation(data)

#' #### Diagnostics
res.assagwt_w %>% DIC() # 5882.934
res.assagwt_w %>% testSSEF()
res.assagwt_w %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagwt_w %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagwt_w %>% plotAC(var.regex)

#' #### Model
res.assagwt_w$model

#' #### Results
res.assagwt_w



#' ### res.assaghwtw {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + hospital + wardtype + ward
#'

# run_model(data, "assaghwtw", directory, thin = 20)
res.assaghwtw <- get_model("assaghwtw", directory)

#' #### Posterior
#+ res.assaghwtw, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assaghwtw %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assaghwtw %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assaghwtw %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assaghwtw %>% plotcorrelation(data)

#' #### Diagnostics
res.assaghwtw %>% DIC() # 5885.554
res.assaghwtw %>% testSSEF()
res.assaghwtw %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assaghwtw %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assaghwtw %>% plotAC(var.regex)

#' #### Model
res.assaghwtw$model

#' #### Results
res.assaghwtw




#' ### res.assagh_wt_w {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + hospital_{wardtype_{ward}}
#'

# run_model(data, "assagh_wt_w", directory, thin = 20)
res.assagh_wt_w <- get_model("assagh_wt_w", directory)

#' #### Posterior
#+ res.assagh_wt_w, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagh_wt_w %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagh_wt_w %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagh_wt_w %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagh_wt_w %>% plotcorrelation(data)

#' #### Diagnostics
res.assagh_wt_w %>% DIC() # 5882.687
res.assagh_wt_w %>% testSSEF()
res.assagh_wt_w %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagh_wt_w %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagh_wt_w %>% plotAC(var.regex)

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
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward + clinical
#'

# run_model(data, "assagwc", directory, thin = 20)
res.assagwc <- get_model("assagwc", directory)

#' #### Posterior
#+ res.assagwc, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagwc %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagwc %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwc %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwc %>% plotcorrelation(data)

#' #### Diagnostics
res.assagwc %>% DIC() # 5886.387
res.assagwc %>% testSSEF()
res.assagwc %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagwc %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagwc %>% plotAC(var.regex)

#' #### Model
res.assagwc$model

#' #### Results
res.assagwc




#' ### res.assagwcst {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward + clinical + sample_type
#'

# run_model(data, "assagwcst", directory, thin = 30)
res.assagwcst <- get_model("assagwcst", directory)

#' #### Posterior
#+ res.assagwcst, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagwcst %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagwcst %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwcst %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwcst %>% plotcorrelation(data)

#' #### Diagnostics
res.assagwcst %>% DIC() # 5835.189
res.assagwcst %>% testSSEF()
res.assagwcst %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagwcst %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagwcst %>% plotAC(var.regex)

#' #### Model
res.assagwcst$model

#' #### Results
res.assagwcst




#' ### res.assagwc_st {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward + clinical_{sampletype}
#'

# run_model(data, "assagwc_st", directory, thin = 30)
res.assagwc_st <- get_model("assagwc_st", directory)

#' #### Posterior
#+ res.assagwc_st, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.assagwc_st %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.assagwc_st %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.assagwc_st %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.assagwc_st %>% plotcorrelation(data)

#' #### Diagnostics
res.assagwc_st %>% DIC() # 5846.981
res.assagwc_st %>% testSSEF()
res.assagwc_st %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.assagwc_st %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.assagwc_st %>% plotAC(var.regex)

#' #### Model
res.assagwc_st$model

#' #### Results
res.assagwc_st



#' ### res.a_cssagwst {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup + ward + sampletype
#'

# run_model(data, "a_cssagwst", directory, thin = 10)
res.a_cssagwst <- get_model("a_cssagwst", directory)

#' #### Posterior
#+ res.a_cssagwst, fig.height = 10
var.regex <- "(a.prob)|(prob.of)|(intercept)|(sd)"
params <- get_params()
labels <- get_labels(data)
res.a_cssagwst %>% densityplot(data, var.regex, params, labels)

#' #### Statistics
res.a_cssagwst %>% summarise_samples(data)

#' #### AMR Summary
#+ fig.height = 10
res.a_cssagwst %>% antibioticsplot(data)

#' #### AMR Correlation
#+ fig.height = 10
res.a_cssagwst %>% plotcorrelation(data)

#' #### Diagnostics
res.a_cssagwst %>% DIC() # 5773.768
res.a_cssagwst %>% testSSEF()
res.a_cssagwst %>% testPSRF()

#' #### Traceplot
#+ fig.height = 6
res.a_cssagwst %>% traceplot(var.regex)

#' #### Autocorrelation
#+ fig.height = 6
res.a_cssagwst %>% plotAC(var.regex)

#' #### Model
res.a_cssagwst$model

#' #### Results
res.a_cssagwst


#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))
DICtable(c("res.assg", "res.assag", "res.assag2", "res.assage",
           "res.assagesq", "res.assagg"))
DICtable(c("res.assagh", "res.assagwt", "res.assagw", "res.assagwt_w",
           "res.assaghwtw", "res.assagh_wt_w"))
DICtable(c("res.assagwc", "res.assagwcst", "res.assagwc_st", "res.a_cssagwst"))










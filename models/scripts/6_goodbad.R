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
#'     toc_collapsed: true
#' ---
#' Last updated: `r Sys.time()`

#+ setup, include=F

# Split posterior plots into good group carriage and good group clinical etc.
# (include indices)
# Add text to report regarding differences between ac models and age groupings
#

library(SpARKjags)
library(runjags)
library(dplyr)
set.seed(1234)

knitr::opts_chunk$set(warning = FALSE)

# Find models here (installed with the package)
directory <- "goodbad_models"
# Save results here (use absolute path so that results can be found when
# RMarkdown report is generated)
res_dir <- file.path("", "Users", "Soniam", "Desktop", "git", "SpARK",
                     "SpARKjags", "results", directory)
model_dir <- file.path("", "Users", "Soniam", "Desktop", "git", "SpARK",
                       "SpARKjags", "models", directory)

data_full <- jags_data(classification = "all",
                       categories = "human",
                       pathogen = "Klebsiella pneumoniae",
                       removeQuinPen = F)

# Remove Quinolone and Penicillin as levels in antibiotic_classes
data <- jags_data(classification = "all",
                  categories = "human",
                  pathogen = "Klebsiella pneumoniae",
                  removeQuinPen = T)

# data_ah <- jags_data(classification = "all",
#                      categories = c("human", "animal"),
#                      pathogen = "Klebsiella pneumoniae",
#                      removeQuinPen = T)

#+


#' ### 0. Data summary
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

# a_naive -----------------------------------------------------------------

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
#' ## 1.1. Experimental structure
#'
#' ### a_naive {.tabset}
#' response ~ antibiotic_class (naive model)
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(model_dir, "a_naive.R"),
                            save_to = res_dir)
res.a_naive <- get_model(path)
model_name <- "res.a_naive"

#' #### Posterior
#+ res.a_naive, fig.height = 6
plot_density(model = get(model_name),
             data = data,
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Diagnostics
DIC(model_name = get(model_name))
testSSEF(model = get(model_name))
testPSRF(model = get(model_name))

#' #### Trace plot
#+ fig.height = 6
plot_caterpillar(model = get(model_name),
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)

# a -----------------------------------------------------------------------

#' ### a {.tabset}
#' $$y_{ij} ~ Bern(a_{goodbad})$$
#' where $a$ is the antibiotic class, $j$

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(model_dir, "a.R"),
                            save_to = res_dir,
                            thin = 10)
res.a <- get_model(path)
model_name <- "res.a"

#' #### Posterior
#+ res.a, fig.height = 8
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)

# a_minus -----------------------------------------------------------------

#' ### a_minus {.tabset}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(model_dir, "a_minus.R"),
                            save_to = res_dir)
res.a_minus <- get_model(path)
model_name <- "res.a_minus"

#' #### Posterior
#+ res.a_minus, fig.height = 8
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)

# a_plus ------------------------------------------------------------------

#' ### a_plus {.tabset}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(model_dir, "a_plus.R"),
                            save_to = res_dir)
res.a_plus <- get_model(path)
model_name <- "res.a_plus"

#' #### Posterior
#+ res.a_plus, fig.height = 8
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)

# ac1 ---------------------------------------------------------------------

#' ### ac1 {.tabset}
#' Combination of ac1
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(model_dir, "ac1.R"),
                            save_to = res_dir)
res.ac1 <- get_model(path)
model_name <- "res.ac1"

#' #### Posterior
#+ res.ac1, fig.height = 8
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)

#' ### Summary of DIC results {-}
#'
DICtable(c("res.a_naive", "res.a", "res.a_minus", "res.a_plus", "res.ac1"))

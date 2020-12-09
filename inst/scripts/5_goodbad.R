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
                       "SpARKjags", "inst", directory)

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

# Louise's test data
test_data <- get_test_data()

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
                            SpARKjags_model = file.path(directory, "a_naive.R"),
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




#' ### a {.tabset}
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

data_ah$data.animal$response %>% data.frame() %>%
  select(Carbapenem) %>% sum(na.rm = T)



# -------------------------------------------------------------------------

#' # Louise's test data

#' ### test_ac1 {.tabset}
#' response ~ antibiotic.class_{goodbad}\
#' goodbad ~ clinical
#'
#' For a particular sample, the probability of resistance (response variable)
#' to a particular antibiotic class is dependent on whether the sample is
#' resistance to other antibiotic classes, as well as the probability of that
#' sample belonging to the good or the bad group, where samples are treated
#' independently from each location (hospital, gp, outpatient, or volunteer).
#'
#' The probability of belonging to the good or the bad group is dependent on
#' whether the sample is clinical or carriage.
#'

path <- run_SpARKjags_model(data = test_data,
                            SpARKjags_model = file.path(directory, "test_ac1.R"),
                            save_to = res_dir)
res.test_ac1 <- get_model(path)
model_name <- "res.test_ac1"

#' #### Posterior
#+ res.test_ac1, fig.height = 10
plot_density(model = get(model_name),
             data = test_data,
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = test_data)

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




#' ### test_ac2 {.tabset}
#'

path <- run_SpARKjags_model(data = test_data,
                            SpARKjags_model = file.path(directory, "test_ac2.R"),
                            save_to = res_dir)
res.test_ac2 <- get_model(path)
model_name <- "res.test_ac2"

#' #### Posterior
#+ res.test_ac2, fig.height = 10
plot_density2(model = get(model_name),
              data = test_data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
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




#' ### test_ac3 {.tabset}
#'

path <- run_SpARKjags_model(data = test_data,
                            SpARKjags_model = file.path(directory, "test_ac3.R"),
                            save_to = res_dir)
res.test_ac3 <- get_model(path)
model_name <- "res.test_ac3"

#' #### Posterior
#+ res.test_ac3, fig.height = 10
plot_density2(model = get(model_name),
              data = test_data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
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




#' ### test_ac4 {.tabset}
#'

path <- run_SpARKjags_model(data = test_data,
                            SpARKjags_model = file.path(directory, "test_ac4.R"),
                            save_to = res_dir)
res.test_ac4 <- get_model(path)
model_name <- "res.test_ac4"

#' #### Posterior
#+ res.test_ac4, fig.height = 10
plot_density2(model = get(model_name),
              data = test_data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
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




#' ### test_a_c {.tabset}
#'

path <- run_SpARKjags_model(data = test_data,
                            SpARKjags_model = file.path(directory, "test_a_c.R"),
                            save_to = res_dir)
res.test_a_c <- get_model(path)
model_name <- "res.test_a_c"

#' #### Posterior
#+ res.test_a_c, fig.height = 10
plot_density2(model = get(model_name),
              data = test_data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
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


DICtable(c("res.test_ac1", "res.test_ac2", "res.test_ac3", "res.test_ac4",
           "res.test_a_c"))




# -------------------------------------------------------------------------

#' ### ac1 {.tabset}
#' response ~ antibiotic.class_{goodbad}\
#' goodbad ~ clinical
#'
#' For a particular sample, the probability of resistance (response variable)
#' to a particular antibiotic class is dependent on whether the sample is
#' resistance to other antibiotic classes, as well as the probability of that
#' sample belonging to the good or the bad group, where samples are treated
#' independently from each location (hospital, gp, outpatient, or volunteer).
#'
#' The probability of belonging to the good or the bad group is dependent on
#' whether the sample is clinical or carriage.
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




#' ### ac2 {.tabset}
#' response ~ antibiotic.class_{goodbad} + clinical
#'
#' For a particular sample, the probability of resistance (response variable)
#' to a particular antibiotic class is dependent on whether the sample is
#' resistance to other antibiotic classes, as well as the probability of that
#' sample belonging to the good or the bad group, and whether the sample is
#' clinical or carriage, where samples are treated independently from each
#' location (hospital, gp, outpatient, or volunteer).
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ac2.R"),
                            save_to = res_dir)
res.ac2 <- get_model(path)
model_name <- "res.ac2"

#' #### Posterior
#+ res.ac2, fig.height = 10
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




#' ### ac3 {.tabset}
#' Combination of ac1 and ac2
#' response ~ antibiotic.class_{goodbad} + clinical\
#' goodbad ~ clinical
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ac3.R"),
                            save_to = res_dir)
res.ac3 <- get_model(path)
model_name <- "res.ac3"

#' #### Posterior
#+ res.ac3, fig.height = 10
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




#' ### ac4 {.tabset}
#' Combination of ac1 and a_c
#' response ~ antibiotic.class_{goodbad,clinical}\
#' goodbad ~ clinical\
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "ac4.R"),
                            save_to = res_dir)
res.ac4 <- get_model(path)
model_name <- "res.ac4"

#' #### Posterior
#+ res.ac4, fig.height = 10
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




#' ### a_c {.tabset}
#' response ~ antibiotic.class_{goodbad,clinical}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "a_c.R"),
                            save_to = res_dir)
res.a_c <- get_model(path)
model_name <- "res.a_c"

#' #### Posterior
#+ res.a_c, fig.height = 10
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
DICtable(c("res.a_naive", "res.a", "res.ac1", "res.ac2", "res.a_c","res.ac3",
           "res.ac4"))




#' ### (a_c)sm {.tabset}
#' response ~ antibiotic.class_{goodbad,clinical} + sample.month
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)sm.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_csm <- get_model(path)
model_name <- "res.a_csm"
file_name <- "res.(a_c)sm"

#' #### Posterior
#+ res.a_csm, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ss.R"),
                            save_to = res_dir,
                            thin = 30)
res.a_css <- get_model(path)
file_name <- "res.(a_c)ss"
model_name <- "res.a_css"

#' #### Posterior
#+ res.a_css, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)


#' ### Summary of DIC results {-}
#'

DICtable(c("res.a_naive", "res.a", "res.ac1", "res.ac2", "res.a_c","res.ac3",
           "res.ac4"))
DICtable(c("res.a_csm", "res.a_css"))


#'
#' ## 1.2. Demographic structure
#' Find the best model with gender and age
#'
#' ### (a_c)ssg {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + gender
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ssg.R"),
                            save_to = res_dir,
                            thin = 30)
res.a_cssg <- get_model(path)
file_name <- "res.(a_c)ssg"
model_name <- "res.a_cssg"

#' #### Posterior
#+ res.a_cssg, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ssag {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ssag.R"),
                            save_to = res_dir,
                            thin = 30)
res.a_cssag <- get_model(path)
file_name <- "res.(a_c)ssag"
model_name <- "res.a_cssag"

#' #### Posterior
#+ res.a_cssag, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ssag2 {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup2
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ssag2.R"),
                            save_to = res_dir,
                            thin = 30)
res.a_cssag2 <- get_model(path)
file_name <- "res.(a_c)ssag2"
model_name <- "res.a_cssag2"

#' #### Posterior
#+ res.a_cssag2, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ssage {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sampling_month + age
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ssage.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssage <- get_model(path)
file_name <- "res.(a_c)ssage"
model_name <- "res.a_cssage"

#' #### Posterior
#+ res.a_cssage, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ssagesq {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + age^2 + age
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ssagesq.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssagesq <- get_model(path)
file_name <- "res.(a_c)ssagesq"
model_name <- "res.a_cssagesq"

#' #### Posterior
#+ res.a_cssagesq, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)g {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup + gender
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "(a_c)ss(ag2)g.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2g <- get_model(path)
file_name <- "res.(a_c)ss(ag2)g"
model_name <- "res.a_cssag2g"

#' #### Posterior
#+ res.a_cssag2g, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)\_g {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup_gender
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "(a_c)ss(ag2)_g.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2_g <- get_model(path)
file_name <- "res.(a_c)ss(ag2)_g"
model_name <- "res.a_cssag2_g"

#' #### Posterior
#+ res.a_cssag2_g, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)



#' ### Summary of DIC results {-}
#'

DICtable(c("res.a_naive", "res.a", "res.ac1", "res.ac2", "res.a_c","res.ac3",
           "res.ac4"))
DICtable(c("res.a_csm", "res.a_css"))
DICtable(c("res.a_cssg", "res.a_cssag", "res.a_cssag2", "res.a_cssage",
           "res.a_cssagesq", "res.a_cssag2g", "res.a_cssag2_g"))

#'
#' ## 1.3. Hospital structure
#' Find the best model with hospital, wardtype, and ward
#'
#' ### (a_c)ss(ag2)h {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup + hospital
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ss(ag2)h.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2h <- get_model(path)
file_name <- "res.(a_c)ss(ag2)h"
model_name <- "res.a_cssag2h"

#' #### Posterior
#+ res.a_cssag2h, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)wt {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup + wardtype
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "(a_c)ss(ag2)wt.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2wt <- get_model(path)
file_name <- "res.(a_c)ss(ag2)wt"
model_name <- "res.a_cssag2wt"

#' #### Posterior
#+ res.a_cssag2wt, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)w {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup + ward
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory, "(a_c)ss(ag2)w.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2w <- get_model(path)
file_name <- "res.(a_c)ss(ag2)w"
model_name <- "res.a_cssag2w"

#' #### Posterior
#+ res.a_cssag2w, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)wt_w {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup +
#' wardtype_{ward}
#'

path <- run_SpARKjags_model(data = data,
                            SpARKjags_model = file.path(directory,
                                                        "(a_c)ss(ag2)wt_w.R"),
                            save_to = res_dir,
                            thin = 20)
res.a_cssag2wt_w <- get_model(path)
file_name <- "res.(a_c)ss(ag2)wt_w"
model_name <- "res.a_cssag2wt_w"

#' #### Posterior
#+ res.a_cssag2wt_w, fig.height = 10
plot_density2(model = get(model_name),
              data = data,
              save_to = file.path(res_dir, "density_plots"),
              filename = paste0(gsub("res.", "", file_name), ".rds"))

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
                 filename = paste0(gsub("res.", "", file_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
#+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", file_name), ".R"))), eval = FALSE, class.source = 'fold-show'

#' #### Results
get(model_name)




#' ### (a_c)ss(ag2)h_wt_w {.tabset}
#' response ~ antibiotic_class_{goodbad,clinical} + sample_season + agegroup +
#' hospital_{wardtype_{ward}}
#'

#' path <- run_SpARKjags_model(data = data,
#'                             SpARKjags_model = file.path(directory,
#'                                                         "(a_c)ss(ag2)h_wt_w.R"),
#'                             save_to = res_dir,
#'                             thin = 20)
#' res.a_cssag2h_wt_w <- get_model(path)
#' model_name <- "res.a_cssag2h_wt_w"
#'
#' #' #### Posterior
#' #+ res.a_cssag2h_wt_w, fig.height = 10
#' plot_density2(model = get(model_name),
#'               data = data,
#'               save_to = file.path(res_dir, "density_plots"),
#'               filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Statistics
#' summarise_samples(model = get(model_name),
#'                   data = data)
#'
#' #' #### AMR Summary
#' #+ fig.height = 8, fig.width = 10
#' plot_antibiotics(model = get(model_name),
#'                  data = data)
#'
#' #' #### AMR Correlation
#' #+ fig.height = 10
#' plot_correlation(model = get(model_name),
#'                  data = data)
#'
#' #' #### Diagnostics
#' DIC(model_name = get(model_name))
#' testSSEF(model = get(model_name))
#' testPSRF(model = get(model_name))
#'
#' #' #### Trace plot
#' #+ fig.height = 6
#' plot_caterpillar(model = get(model_name),
#'                  save_to = file.path(res_dir, "caterpillar_plots"),
#'                  filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Autocorrelation
#' #+ fig.height = 6
#' plot_autocorr(model = get(model_name))
#'
#' #' #### Model
#' #+ code = readLines(file.path(model_dir, paste0(gsub("res.", "", model_name), ".R"))), eval = FALSE, class.source = 'fold-show'
#'
#' #' #### Results
#' get(model_name)


#' ### Summary of DIC results {-}
#'

DICtable(c("res.a_naive", "res.a", "res.ac1", "res.ac2", "res.a_c","res.ac3",
           "res.ac4"))
DICtable(c("res.a_csm", "res.a_css"))
DICtable(c("res.a_cssg", "res.a_cssag", "res.a_cssag2", "res.a_cssage",
           "res.a_cssagesq", "res.a_cssag2g", "res.a_cssag2_g"))
DICtable(c("res.a_cssag2h", "res.a_cssag2wt", "res.a_cssag2w",
           "res.a_cssag2wt_w"
           # , "res.a_cssag2h_wt_w"
           ))


#' ### res.assagwcst {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical + sample_type
#'

#' path <- run_SpARKjags_model(data = data,
#'                             SpARKjags_model = file.path(directory, "assagwcst.R"),
#'                             save_to = res_dir,
#'                             thin = 30)
#' res.assagwcst <- get_model(path)
#' model_name <- "res.assagwcst"
#'
#' #' #### Posterior
#' #+ res.assagwcst, fig.height = 10
#' plot_density(model = get(model_name),
#'              data = data,
#'              save_to = file.path(res_dir, "density_plots"),
#'              filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Statistics
#' summarise_samples(model = get(model_name),
#'                   data = data)
#'
#' #' #### AMR Summary
#' #+ fig.height = 8, fig.width = 10
#' plot_antibiotics(model = get(model_name),
#'                  data = data)
#'
#' #' #### AMR Correlation
#' #+ fig.height = 10
#' plot_correlation(model = get(model_name),
#'                  data = data)
#'
#' #' #### Diagnostics
#' DIC(model_name = get(model_name))
#' testSSEF(model = get(model_name))
#' testPSRF(model = get(model_name))
#'
#' #' #### Trace plot
#' #+ fig.height = 6
#' plot_caterpillar(model = get(model_name),
#'                  save_to = file.path(res_dir, "caterpillar_plots"),
#'                  filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Autocorrelation
#' #+ fig.height = 6
#' plot_autocorr(model = get(model_name))
#'
#' #' #### Model
#' get(model_name)$model
#'
#' #' #### Results
#' get(model_name)




#' ### res.assagwc_st {.tabset}
#' response ~ antibiotic_class_{goodbad} + sample_season + agegroup + ward +
#' clinical_{sampletype}
#'

#' path <- run_SpARKjags_model(data = data,
#'                             SpARKjags_model = file.path(directory, "assagwc_st.R"),
#'                             save_to = res_dir,
#'                             thin = 30)
#' res.assagwc_st <- get_model(path)
#' model_name <- "res.assagwc_st"
#'
#' #' #### Posterior
#' #+ res.assagwc_st, fig.height = 10
#' plot_density(model = get(model_name),
#'              data = data,
#'              save_to = file.path(res_dir, "density_plots"),
#'              filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Statistics
#' summarise_samples(model = get(model_name),
#'                   data = data)
#'
#' #' #### AMR Summary
#' #+ fig.height = 8, fig.width = 10
#' plot_antibiotics(model = get(model_name),
#'                  data = data)
#'
#' #' #### AMR Correlation
#' #+ fig.height = 10
#' plot_correlation(model = get(model_name),
#'                  data = data)
#'
#' #' #### Diagnostics
#' DIC(model_name = get(model_name))
#' testSSEF(model = get(model_name))
#' testPSRF(model = get(model_name))
#'
#' #' #### Trace plot
#' #+ fig.height = 6
#' plot_caterpillar(model = get(model_name),
#'                  save_to = file.path(res_dir, "caterpillar_plots"),
#'                  filename = paste0(gsub("res.", "", model_name), ".rds"))
#'
#' #' #### Autocorrelation
#' #+ fig.height = 6
#' plot_autocorr(model = get(model_name))
#'
#' #' #### Model
#' get(model_name)$model
#'
#' #' #### Results
#' get(model_name)




#'
#' ### Summary of DIC results {-}
#'
# DICtable(c("res.a_naive", "res.a", "res.asm", "res.ass"))
# DICtable(c("res.assg", "res.assag", "res.assag2", "res.assage",
#            "res.assagesq", "res.assagg"))
# DICtable(c("res.assagh", "res.assagwt", "res.assagw", "res.assagwt_w",
#            "res.assaghwtw"))
# , "res.assagh_wt_w"))
# DICtable(c("res.assagwc", "res.assagwcst", "res.assagwc_st"))
# "res.a_cssagwst"


#'
#' ## 1.4. Animals
#'

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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




#' ### res.a_wild_subsets {.tabset}
#' response ~ antibiotic.class_{goodbad}
#'

path <- run_SpARKjags_model(data = data_ah,
                            SpARKjags_model = file.path(directory,
                                                        "a_wild_subsets.R"),
                            save_to = res_dir,
                            thin = 10)
res.a_wild_subsets <- get_model(path)
model_name <- "res.a_wild_subsets"

#' #### Posterior
#+ res.a_wild_subsets, fig.height = 10
plot_density(model = get(model_name),
             data = data_ah,
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)




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
             save_to = file.path(res_dir, "density_plots"),
             filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Statistics
summarise_samples(model = get(model_name),
                  data = data_ah)

#' #### AMR Summary
#+ fig.height = 8, fig.width = 10
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
                 save_to = file.path(res_dir, "caterpillar_plots"),
                 filename = paste0(gsub("res.", "", model_name), ".rds"))

#' #### Autocorrelation
#+ fig.height = 6
plot_autocorr(model = get(model_name))

#' #### Model
get(model_name)$model

#' #### Results
get(model_name)



#' ---
#' title: Jags models
#' date: Sep 2019
#' ---
#' Last updated: `r Sys.time()`

#+ r setup, include=FALSE


# output:
#   html_document:
#     code_folding: hide

library(SpARKjags)
library(SpARK)
library(dplyr)
library(magrittr)
library(runjags)
library(ggplot2)

knitr::opts_chunk$set(message = FALSE,
                      warning = F,
                      cache = F,
                      echo = F,
                      results = "hide")

set.seed(1234)



# Generate JAGS data
data.all <- get_data(classification = "Carbapenem",
                     kleb = "Klebsiella pneumoniae")

lookup <- data.all$lookup_tables
counts.all.II <- data.all$counts

data.jags <- list(h_resist = data.all$hospital_dat$class_interpretation,
                  ward = data.all$sample_dat$ward,
                  ward_type = data.all$ward_dat$ward_type,
                  hospital = data.all$ward_dat$hospital,
                  h_sample_GUID = data.all$hospital_dat$sample_GUID,
                  h_bacteria = data.all$hospital_dat$bacteria,
                  gp_resist = data.all$gp_dat$class_interpretation,
                  gp_sample_GUID = data.all$gp_dat$sample_GUID,
                  gp_bacteria = data.all$gp_dat$bacteria,
                  v_resist = data.all$volunteer_dat$class_interpretation,
                  v_sample_GUID = data.all$volunteer_dat$sample_GUID,
                  v_bacteria = data.all$volunteer_dat$bacteria,
                  o_resist = data.all$outpatients_dat$class_interpretation,
                  o_sample_GUID = data.all$outpatients_dat$sample_GUID,
                  o_bacteria = data.all$outpatients_dat$bacteria,
                  clinical = data.all$sample_type_dat$clinical,
                  sample_type = data.all$sample_dat$sample_type,
                  N_patients = data.all$N_hospital_dat,
                  N_hosp = data.all$N_hospital,
                  N_gp = data.all$N_gp_dat,
                  N_volunteers = data.all$N_volunteer_dat,
                  N_outpatients = data.all$N_out_dat,
                  N_sample = data.all$N_sample_dat,
                  N_ward = data.all$N_ward_dat,
                  gender = data.all$sample_dat$gender,
                  age = data.all$sample_dat$age,
                  age_group = data.all$sample_dat$age_group,
                  age_group2 = data.all$sample_dat$age_group2,
                  N_age_group = data.all$N_age_group,
                  N_sample_type = data.all$N_sample_type_dat
)

l_clin <- data.all$lookup_tables$clinical
data.jags$ncarr <- l_clin$index[l_clin$clinical=="no"]
data.jags$nclin <- l_clin$index[l_clin$clinical=="yes"]
data.jags$v_clinical  <- data.jags$ncarr # not
data.jags$gp_clinical <- data.jags$nclin # yes
data.jags$o_clinical  <- data.jags$nclin # yes
l_gen <- data.all$lookup_tables$gender
data.jags$female <- l_gen$index[l_gen$gender=="F"]
data.jags$male <- l_gen$index[l_gen$gender=="M"]
data.jags$genders <- c(data.jags$female, data.jags$male)
data.jags$hosp_wards <- sort(unique(data.jags$ward[data.jags$h_sample_GUID]))
data.jags$hosp_wardtypes <- sort(unique(data.jags$ward_type[data.jags$ward[data.jags$h_sample_GUID]]))
data.jags$agegroups <- sort(unique(data.jags$age_group2))
data.jags$bact_species <-  sort(unique(c(data.jags$h_bacteria,
                                         data.jags$gp_bacteria,
                                         data.jags$v_bacteria,
                                         data.jags$o_bacteria)))


lookup$gender

model <- run.jags("../individual_models/gen1.R",
                  data = data.jags,
                  n.chains = 2)

# Autocorrelation plot for each paramter shows the degree of correlation
# between MCMC samples separated by different lags. For example, a lag of
# 0 represents the degree of correlation between each MCMC sample and itself
# (obviously this will be a correlation of 1). A lag of 1 represents the
# degree of correlation between each MCMC sample and the next sample along
# the Chain and so on. In order to be able to generate unbiased estimates
# of parameters, the MCMC samples should be independent (uncorrelated).
model %>% coda::as.mcmc.list() %>% coda::autocorr.plot()
model %>% coda::as.mcmc.list() %>% coda::crosscorr.plot()


# Gelman and Rubin's convergence diagnostic
model %>% coda::as.mcmc.list() %>% coda::gelman.diag(multivariate = T)
model %>% coda::as.mcmc.list() %>% coda::gelman.plot()

# Get a sample of posterior draws from a model as a tibble
tidybayes::tidy_draws(model)

rjags::coda.samples(as.jags(model), "gender.effect", 10)





# DIC looks at the model comparison based on a single overall estimate of predictive log likelihood. While LOO and WAIC account for subject predictive log-likelohhod across the multiple iterations



# Male samples are have a 14% probability of having carbapenem resistance.
# Female samples are have a 23% probability of having carbapenem resistance.
# Female samples are 1.8 times more likely to have carbapenem resistance than
# male samples.

# get_parameters(results.g)
#
# lab.g <- results.frame(Parameter = paste0("gender.effect[", lookup$gender$index, "]"),
#                        Label = lookup$gender$gender)
# plot_jags(results.g, lab.g$Parameter, lab.g, 0)
# plot_jags(results.g, "sd")


# # Convert log-odds predictions to probabilities
# tmp <- data.frame(results.age$summaries)[, "Mean", drop = F] %>%
#   tibble::rownames_to_column("variable") %>%
#   mutate(odds = case_when(grepl("effect", variable) ~ exp(Mean),
#                           grepl("intercept", variable) ~ exp(Mean)),
#          probability = case_when(grepl("effect", variable) ~ plogis(Mean),
#                                  grepl("intercept", variable) ~ plogis(Mean)))
# tmp
#
# # Predicted resistances from parameter estimates
# intercept <- tmp %>% filter(variable == "intercept") %$% Mean
# beta1 <- tmp %>% filter(variable == "age.effect") %$% Mean
# plot_line <- data.frame(age = 0:100) %>%
#   mutate(y = beta1*age + intercept,
#          prob = plogis(y))
#
# # Extract resistances from dataset
# dots <- rbind.data.frame(
#   data.all$hospital_dat %>% select(GUID, sample_GUID, class_interpretation),
#   data.all$gp_dat %>% select(GUID, sample_GUID, class_interpretation),
#   data.all$outpatients_dat %>% select(GUID, sample_GUID, class_interpretation),
#   data.all$volunteer_dat %>% select(GUID, sample_GUID, class_interpretation)) %>%
#   merge(data.all$sample_dat %>% select(sample_GUID, age)) %>%
#   rename(y = class_interpretation) %>%
#   select(y, age) %>%
#   mutate(breaks = cut(age, breaks = seq(0,100,2), labels = seq(1,99,2),
#                       include.lowest = TRUE),
#          breaks = as.numeric(as.character(breaks)))
# plot_dots <- dots %>%
#   group_by(y, breaks) %>%
#   summarise(n = n()) %>%
#   mutate(pct = ifelse(y == 0, n / sum(n), 1 - n / sum(n)))
#
# # Plot logistic regression curve
# ggplot() + theme_minimal() +
#   geom_segment(data = plot_dots, size = 4, show.legend = FALSE,
#                aes(x = breaks, xend = breaks, y = y, yend = pct,
#                    colour = factor(y))) +
#   geom_rug(aes(x = age), subset(dots, y == 1), sides = "t",
#            colour = "grey50", alpha = 0.4) +
#   geom_rug(aes(x = age), subset(dots, y == 0), sides = "b",
#            colour = "grey50", alpha = 0.4) +
#
#   geom_line(aes(x = age, y = prob), plot_line, colour = "grey50", lwd = 1) +
#   labs(x = "age", y = "y")




























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

library(knitr)
library(runjags)
library(rjags)
library(dplyr)
library(SpARK)
library(SpARKjags)
library(lubridate)
library(ggplot2)
library(broom)

knitr::opts_chunk$set(message = FALSE,
                      warning = F,
                      cache = T,
                      echo = F,
                      results = "hide")




# 1. Generate JAGS data ---------------------------------------------------

# data.all <- get_data(classification = "Carbapenem", kleb = "all")
# data.all <- get_data(classification = "Carbapenem", pathogen = "all")
data.all <- get_data(classification = "Carbapenem", kleb = "Klebsiella pneumoniae")
# data.all <- get_data(classification = "Carbapenem", pathogen = "Klebsiella_pneumoniae")
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
N_bacteria <- length(data.jags$bact_species)

table(data.jags$h_resist,
      data.jags$h_bacteria,
      data.jags$clinical[data.jags$sample_type[data.jags$h_sample_GUID]])

jags.mod.null <- run.jags("Zmodels/null.jags",
                          data = data.jags,
                          n.chains = 2)
#+ results = "markup"
jags.mod.null # 943.9

results.g <- run.jags("models/gender.jags",
                      data = data.jags,
                      n.chains = 2)

results.g # 933.5

results.ag <- run.jags("models/agegroup.jags",
                       data = data.jags,
                       n.chains = 2)
results.ag # 938.8

results.ag2 <- run.jags("models/agegroup2.jags", # --------- !
                        data = data.jags,
                        n.chains = 2)
results.ag2 # 930.9

results.age <- run.jags("models/age.jags", # --------- !
                        data = data.jags,
                        n.chains = 2)
results.age # 943.7

# tmp.age <- model$summaries %>%
#   data.frame() %>%
#   tibble::rownames_to_column() %>%
#   rename(var = rowname,
#          logits_resistance = Mean) %>%
#   mutate(pred_prob_resistance = exp(logits_resistance) /
#            (1 + exp(logits_resistance)))
#
# data.jags

results.age2 <- run.jags("models/age2.jags", # --------- !
                         data = data.jags,
                         n.chains = 2)

results.age2 # 922.0

results.c <- run.jags("models/clinical.jags",
                      data = data.jags,
                      n.chains = 2)
results.c # 899.6

results.h <- run.jags("models/hospital.jags",
                      data = data.jags,
                      n.chains = 2)
results.h # 868.0

results.wt <- run.jags("models/wardtype.jags",
                       data = data.jags,
                       n.chains = 2)
results.wt # 855.2

if (N_bacteria > 1)
{
  results.b <- run.jags("models/bacteria.jags",
                        data = data.jags,
                        n.chains = 2)
  results.b # 1013.8

  results.b_c <- run.jags("models/b+c.jags",
                          data = data.jags,
                          n.chains = 2)
  results.b_c # 1013.8
}

results.hw <- run.jags("models/hw.jags",
                       data = data.jags,
                       n.chains = 2)
results.hw # 830.8

results.w <- run.jags("models/ward.jags",
                      data = data.jags,
                      n.chains = 2)
results.w # 828.8

results.wtw <- run.jags("models/wtw.jags",
                        data = data.jags,
                        n.chains = 2)
results.wtw # 826.0

results.wc <- run.jags("models/wc.jags",
                       data = data.jags,
                       n.chains = 2)
results.wc # 825.5

results.cw <- run.jags("models/cw.jags",
                       data = data.jags,
                       n.chains = 2)
results.cw # 823.0

results.w_c <- run.jags("models/w+c.jags",
                        data = data.jags,
                        n.chains = 2)
results.w_c # 821.3

results.w_c %>% density_plot("odds_c", lookup)

results.cwtw <- run.jags("models/cwtw.jags",
                         data = data.jags,
                         n.chains = 2)
results.cwtw # 819.4

results.cwtw %>%
  caterpillar(c("c_diff", "sd.wt", "sd.w"))

results.cwtw %>%
  coda::as.mcmc.list() %>%
  ggmcmc::ggs() %$%
  Parameter %>%
  unique %>% as.character

results.cwtw %>%
  density_plot(c("c_diff", "sd.wt", "sd.w"), lookup) +
  ggplot2::labs(x = "Values")

results.cwtw %>%
  density_plot(c("gp_prob", "o_prob", "v_prob"), lookup)

results.cwtw %>%
  caterpillar("c_diff")

results.cwtwg <- run.jags("models/cwtwg.jags",
                          data = data.jags,
                          n.chains = 2)
results.cwtwg # 810.0

if (N_bacteria > 1)
{
  results.cwtwb <- run.jags("models/cwtwb.jags",
                            data = data.jags,
                            n.chains = 2)
  results.cwtwb # 863.5
}

results.cwtwga <- run.jags("models/cwtwga.jags",
                           data = data.jags,
                           n.chains = 2)
results.cwtwga # 810.0

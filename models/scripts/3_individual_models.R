#' ---
#' title: Seemingly unrelated regressions
#' date: Sep 2019
#' ---
#' Last updated: `r Sys.time()`

#' Here we're modelling resistance for each antibiotic class.
#' When the response variable is carbapenem resistance, use the whole dataset.
#' When it's a different antibiotic class, remove all samples that have
#' carbapenem resistance.

#+ r setup, include=FALSE

library(SpARKcarbapenem)
library(SpARK)
library(dplyr)
library(runjags)
library(ggplot2)
library(coda)

knitr::opts_chunk$set(message = FALSE,
                      warning = F,
                      cache = F,
                      echo = F,
                      results = "hide")

set.seed(1234)

antibiotic_classes <- c("Aminoglycoside", "Penicillin Combination",
                        # "Penicillin",
                        "Monobactam", "Cephalosporin", "Fluoroquinolone",
                        "Colistin", "Carbapenem", "Fosfomycin",
                        "Penicillin (Penams)",
                        # "Quinolone",
                        "Nitrofurantoin",
                        "Tetracycline", "Trimethoprim",
                        "Trimethoprim/Sulfamethoxazole")

for(x in seq_along(antibiotic_classes)) {
  cat("\n\n\nAntibiotic class", x, "of", length(antibiotic_classes), "...")

  if(antibiotic_classes[x] == "Carbapenem") {
    data.jags <- jags_data(classification = antibiotic_classes[x],
                           categories = "human",
                           kleb = "Klebsiella pneumoniae")$jags
  } else {
    data.jags <- jags_data(classification = antibiotic_classes[x],
                           categories = "human",
                           kleb = "Klebsiella pneumoniae",
                           removeCarbapenem = T)$jags
  }

  # If all samples are resistant (or all are not resistant), skip
  if(length(unique(c(unique(data.jags$h_resist),
                     unique(data.jags$gp_resist),
                     unique(data.jags$v_resist),
                     unique(data.jags$o_resist)))) == 1) next

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept
  # intercept ~ N(0, tau)
  # tau       ~ G(0.001, 0.001)
  results.null <- run.jags("inst/individual_models/null.R",
                           data = data.jags,
                           n.chains = 2) # 973.0

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = gender.effect
  # gender.effect ~ N(0, tau)
  # tau           ~ G(0.001, 0.001)
  results.g <- run.jags("inst/individual_models/g.R",
                        data = data.jags,
                        n.chains = 2) # 962.9

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + agegroup.effect
  # agegroup.effect ~ N(0, tau)
  # intercept       ~ N(0, 0.0001)
  # tau             ~ G(0.001, 0.001)
  results.agegroup <- run.jags("inst/individual_models/ag.R",
                               data = data.jags,
                               n.chains = 2) # 966.8

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + agegroup2.effect
  # agegroup2.effect ~ N(0, tau)
  # intercept       ~ N(0, 0.0001)
  # tau             ~ G(0.001, 0.001)
  results.agegroup2 <- run.jags("inst/individual_models/ag2.R",
                                data = data.jags,
                                n.chains = 2) # 962.7 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.effect * AGE_i
  # age.effect ~ N(0, 0.001)
  # intercept  ~ N(0, 0.001)
  results.age <- run.jags("inst/individual_models/age.R",
                          data = data.jags,
                          n.chains = 2) # 973.1 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.sq.effect * AGE^2_i + age.effect * AGE_i
  # age.sq.effect ~ N(0, 0.001)
  # age.effect    ~ N(0, 0.001)
  # intercept     ~ N(0, 0.001)
  results.agesq <- run.jags("inst/individual_models/agesq.R",
                            data = data.jags,
                            n.chains = 2) # 949.0 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.sq.effect * AGE^2_i + age.effect * AGE_i + gender.effect
  # age.effect    ~ N(0, 0.001)
  # age.sq.effect ~ N(0, 0.001)
  # gender.effect ~ N(0, 0.001)
  # intercept     ~ N(0, 0.001)
  results.agesqg <- run.jags("inst/individual_models/agesqg.R",
                             data = data.jags,
                             n.chains = 2) # 939.5 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.sq.effect_{gender} * age^2_i +
  #         age.effect_{gender} * age_i
  # age.effect_{gender}    ~ N(0, 0.001)
  # age.sq.effect_{gender} ~ N(0, 0.001)
  # intercept              ~ N(0, 0.001)
  results.agesq_g <- run.jags("inst/individual_models/agesq_g.R",
                              data = data.jags,
                              n.chains = 2) #  (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.sq.effect_{clin} * age^2_i + age.effect_{clin} * age_i
  # age.effect_{clin}    ~ N(0, 0.001)
  # age.sq.effect_{clin} ~ N(0, 0.001)
  # intercept            ~ N(0, 0.001)
  results.agesq_c <- run.jags("inst/individual_models/agesq_c.R",
                              data = data.jags,
                              n.chains = 2) #  (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + age.sq.effect_{clin} * age^2_i + age.effect_{clin} * age_i
  # age.effect_{clin}    ~ N(0, 0.001)
  # age.sq.effect_{clin} ~ N(0, 0.001)
  # intercept            ~ N(0, 0.001)
  results.agesq_wt <- run.jags("inst/individual_models/agesq_wt.R",
                               data = data.jags,
                               n.chains = 2) #  (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = clinical.effect
  # clinical.effect ~ N(0, tau)
  # tau             ~ G(0.001, 0.001)
  results.c <- run.jags("inst/individual_models/c.R",
                        data = data.jags,
                        n.chains = 2) # 925.2

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,i})
  # p_{sampletype,i} = logit^{-1}(y_{sampletype,i})
  # y_{sampletype,i} = clinical.effect_{sampletype}
  # clinical.effect_{sampletype} ~ N(clinical.effect, tau.stc)
  # clinical.effect              ~ N(0, tau.c)
  # tau.c                        ~ G(0.001, 0.001)
  # tau.stc                      ~ G(0.001, 0.001)
  results.c_st <- run.jags("inst/individual_models/c_st.R",
                           data = data.jags,
                           n.chains = 2) # 888.0 (sonia)


  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clinical.effect + sampletype.effect
  # clinical.effect   ~ N(0, 0.001)
  # sampletype.effect ~ N(0, tau)
  # intercept         ~ N(0, 0.001)
  # tau               ~ G(0.001, 0.001)
  results.cst <- run.jags("inst/individual_models/cst.R",
                          data = data.jags,
                          n.chains = 2) # 887.6 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + sampletype.effect
  # sampletype.effect ~ N(0, tau)
  # intercept         ~ N(0, 0.001)
  # tau               ~ G(0.001, 0.001)
  results.st <- run.jags("inst/individual_models/st.R",
                         data = data.jags,
                         n.chains = 2) # 886.6 (sonia)


  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + hospital.effect
  # hospital.effect ~ N(0, tau)
  # intercept       ~ N(0, 0.001)
  # tau             ~ G(0.001, 0.001)
  results.h <- run.jags("inst/individual_models/h.R",
                        data = data.jags,
                        n.chains = 2) # 893.8

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + wardtype.effect
  # wardtype.effect ~ N(0, tau)
  # intercept       ~ N(0, 0.001)
  # tau             ~ G(0.001, 0.001)
  results.wt <- run.jags("inst/individual_models/wt.R",
                         data = data.jags,
                         n.chains = 2) # 879.9

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = intercept + hosp.effect_{ward}
  # hosp.effect_{ward} ~ N(hosp.effect, tau.ward)
  # hosp.effect        ~ N(0, tau.hosp)
  # intercept          ~ N(0, 0.001)
  # tau.hosp           ~ G(0.001, 0.001)
  # tau.ward           ~ G(0.001, 0.001)
  results.h_w <- run.jags("inst/individual_models/h_w.R",
                          data = data.jags,
                          n.chains = 2) # 852.1

  # RESISTANCE_i ~ Bernoulli(p_{wardtype,i})
  # p_{wardtype,i} = logit^{-1}(y_{wardtype,i})
  # y_{wardtype,i} = intercept + hosp.effect_{wardtype}
  # hosp.effect_{wardtype} ~ N(hosp.effect, tau.wt)
  # hosp.effect        ~ N(0, tau.hosp)
  # intercept          ~ N(0, 0.001)
  # tau.hosp           ~ G(0.001, 0.001)
  # tau.wt             ~ G(0.001, 0.001)
  results.h_wt <- run.jags("inst/individual_models/h_wt.R",
                           data = data.jags,
                           n.chains = 2) # 879.3 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,i})
  # p_{sampletype,i} = logit^{-1}(y_{sampletype,i})
  # y_{sampletype,i} = intercept + hosp.effect_{sampletype}
  # hosp.effect_{sampletype} ~ N(hosp.effect, tau.wt)
  # hosp.effect        ~ N(0, tau.hosp)
  # intercept          ~ N(0, 0.001)
  # tau.hosp           ~ G(0.001, 0.001)
  # tau.st             ~ G(0.001, 0.001)
  results.h_st <- run.jags("inst/individual_models/h_st.R",
                           data = data.jags,
                           n.chains = 2) # 863.7 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + ward.effect
  # ward.effect ~ N(0, tau)
  # intercept   ~ N(0, 0.001)
  # tau         ~ G(0.001, 0.001)
  results.w <- run.jags("inst/individual_models/ward.R",
                        data = data.jags,
                        n.chains = 2) # 850.0

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = intercept + wardtype.effect_{ward}
  # wardtype.effect_{ward} ~ N(wardtype.effect, tau.ward)
  # wardtype.effect        ~ N(0, tau.wt)
  # intercept              ~ N(0, 0.001)
  # tau.wt                 ~ G(0.001, 0.001)
  # tau.ward               ~ G(0.001, 0.001)
  results.wt_w <- run.jags("inst/individual_models/wt_w.R",
                           data = data.jags,
                           n.chains = 2) # 847.7

  # RESISTANCE_i ~ Bernoulli(p_{clinical,i})
  # p_{clinical,i} = logit^{-1}(y_{clinical,i})
  # y_{clinical,i} = intercept + ward.effect_{clinical}
  # ward.effect_{clinical} ~ N(ward.effect, tau.wc)
  # ward.effect            ~ N(intercept, tau.ward)
  # intercept              ~ N(0, 0.001)
  # tau.ward               ~ G(0.001, 0.001)
  # tau.wc                 ~ G(0.001, 0.001)
  results.w_c <- run.jags("inst/individual_models/w_c.R",
                          data = data.jags,
                          n.chains = 2) # 845.2

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,i})
  # p_{sampletype,i} = logit^{-1}(y_{sampletype,i})
  # y_{sampletype,i} = intercept + ward.effect_{sampletype}
  # ward.effect_{sampletype} ~ N(ward.effect, tau.wc)
  # ward.effect            ~ N(intercept, tau.ward)
  # intercept              ~ N(0, 0.001)
  # tau.ward               ~ G(0.001, 0.001)
  # tau.wst                ~ G(0.001, 0.001)
  results.w_st <- run.jags("inst/individual_models/w_st.R",
                           data = data.jags,
                           n.chains = 2) # 844.4 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = clinical.effect_{ward}
  # clinical.effect_{ward} ~ N(clinical.effect, tau)
  # clinical.effect        ~ N(0, 0.001)
  # tau                    ~ G(0.001, 0.001)
  results.c_w <- run.jags("inst/individual_models/c_w.R",
                          data = data.jags,
                          n.chains = 2) # 842.3
  # !!! should I put a nh.effect in here with precision tau.clin? ----
  # !!! no intercept.. is that ok?

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clinical.effect + ward.effect
  # clinical.effect ~ N(0, 0.001)
  # ward.effect     ~ N(0, tau)
  # intercept       ~ N(0, 0.001)
  # tau             ~ G(0.001, 0.001)
  results.wc <- run.jags("inst/individual_models/wc.R",
                         data = data.jags,
                         n.chains = 2) # 840.8
  # should ward.effect have the same tau as gp.effect? ----
  # check c+wt+w too...

  # RESISTANCE_i ~ Bernoulli(p_{wardtype,i})
  # p_{wardtype,i} = logit^{-1}(y_{wardtype,i})
  # y_{wardtype,i} = clinical.effect_{wardtype}
  # clinical.effect_{wardtype} ~ N(clinical.effect, tau)
  # clinical.effect            ~ N(0, 0.001)
  # tau                        ~ G(0.001, 0.001)
  results.c_wt <- run.jags("inst/individual_models/c_wt.R",
                           data = data.jags,
                           n.chains = 2) # 863.5 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = clinical.effect_{wardtype,ward}
  # clinical.effect_{wardtype,ward} ~ N(clinical.effect_{wardtype}, tau.w)
  # clinical.effect_{wardtype}      ~ N(clinical.effect, tau.wt)
  # clinical.effect                 ~ N(0, 0.001)
  # tau.wt                          ~ G(0.001, 0.001)
  # tau.w                           ~ G(0.001, 0.001)
  results.c_wt_w <- run.jags("inst/individual_models/c_wt_w.R",
                             data = data.jags,
                             n.chains = 2) # 842.5 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,i})
  # p_{sampletype,i} = logit^{-1}(y_{sampletype,i})
  # y_{sampletype,i} = clinical.effect_{wardtype,sampletype}
  # clinical.effect_{wardtype,sampletype} ~ N(clinical.effect_{wardtype}, tau.st)
  # clinical.effect_{wardtype}            ~ N(clinical.effect, tau.wt)
  # clinical.effect                       ~ N(0, 0.001)
  # tau.st                                ~ G(0.001, 0.001)
  # tau.wt                                ~ G(0.001, 0.001)
  results.c_wt_st <- run.jags("inst/individual_models/c_wt_st.R",
                              data = data.jags,
                              n.chains = 2) # 848.14 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,i})
  # p_{sampletype,i} = logit^{-1}(y_{sampletype,i})
  # y_{sampletype,i} = clinical.effect_{wardtype,sampletype,ward}
  # clinical.effect_{wardtype,sampletype,ward} ~
  #                             N(clinical.effect_{wardtype,sampletype}, tau.w)
  # clinical.effect_{wardtype,sampletype} ~ N(clinical.effect_{wardtype}, tau.st)
  # clinical.effect_{wardtype}            ~ N(clinical.effect, tau.wt)
  # clinical.effect                       ~ N(0, 0.001)
  # tau.st                                ~ G(0.001, 0.001)
  # tau.wt                                ~ G(0.001, 0.001)
  results.c_wt_st_w <- run.jags("inst/individual_models/c_wt_st_w.R",
                                data = data.jags,
                                n.chains = 2) # 839.75 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = intercept + clinical.effect + wardtype.effect + ward.effect
  # clinical.effect ~ N(0, 0.001)
  # wardtype.effect ~ N(intercept, tau.wt)
  # ward.effect     ~ N(0, tau.w)
  # intercept       ~ N(0, 0.001)
  # tau.wt          ~ G(0.001, 0.001)
  # tau.w           ~ G(0.001, 0.001)
  results.cwtw <- run.jags("inst/individual_models/cwtw.R",
                           data = data.jags,
                           n.chains = 2) # 839.1 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = wtc.effect_{ward,i}        # should i put intercept here? ----
  # wtc.effect_{ward} ~ N(wtc.effect, tau.w)
  # wtc.effect = intercept + clinical.effect + wardtype.effect
  # clinical.effect ~ N(0, 0.001)
  # wardtype.effect ~ N(0, tau.wt)
  # intercept       ~ N(0, 0.001)
  # tau.wt          ~ G(0.001, 0.001)
  # tau.w           ~ G(0.001, 0.001)
  results.cwt_w <- run.jags("inst/individual_models/cwt_w.R",
                            data = data.jags,
                            n.chains = 2) # 837.7

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clinical.effect + wardtype.effect_{ward}
  # wardtype.effect_{ward} ~ N(wardtype.effect, tau.ward)
  # wardtype.effect ~ N(intercept, tau.wt)
  # intercept ~ N(0, 0.001)
  # tau.wt ~ G(0.001, 0.001)
  # tau.w ~ G(0.001, 0.001)
  results.c.wt_w <- run.jags("inst/individual_models/c-wt_w.R",
                             data = data.jags,
                             n.chains = 2) # 836.8 (sonia)

  # I made some edits to the model above !!! ----
  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = wtc.effect_{ward,i}
  # wtc.effect_{ward} ~ N(wtc.effect, tau.w)
  # wtc.effect = intercept + clinical.effect + wardtype.effect
  # clinical.effect ~ N(0, 0.001)
  # wardtype.effect ~ N(0, tau.wt)
  # intercept       ~ N(0, 0.001)
  # tau.wt          ~ G(0.001, 0.001)
  # tau.w           ~ G(0.001, 0.001)
  results.cwt_w2 <- run.jags("inst/individual_models/cwt_w[2].R",
                             data = data.jags,
                             n.chains = 2) # 839.4 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = wtcg.effect_{ward,i}
  # wtcg.effect_{ward} ~ N(wtcg.effect, tau.w)
  # wtcg.effect = intercept + clinical.effect + wardtype.effect + gender.effect
  # clinical.effect ~ N(0, 0.001)
  # gender.effect   ~ N(0, 0.001)
  # wardtype.effect ~ N(0, tau.wt)
  # intercept       ~ N(0, 0.001)
  # tau.wt          ~ G(0.001, 0.001)
  # tau.w           ~ G(0.001, 0.001)
  results.cwtg_w <- run.jags("inst/individual_models/cwtg_w.R",
                             data = data.jags,
                             n.chains = 2) # 830.0

  # RESISTANCE_i ~ Bernoulli(p_{sampletype,ward,i})
  # p_{sampletype,ward,i} = logit^{-1}(y_{sampletype,ward,i})
  # y_{sampletype,ward,i} = intercept + clinical.effect_{sampletype,ward} +
  #                           wardtype.effect_{ward}
  # clinical.effect_{sampletype} ~ N(clinical.effect, tau.cst)
  # clinical.effect              ~ N(0, 0.001)
  # wardtype.effect_{ward}       ~ N(wardtype.effect, tau.w)
  # wardtype.effect              ~ N(0, tau.wt)
  # intercept                    ~ N(0, 0.001)
  # tau.wt                       ~ G(0.001, 0.001)
  # tau.w                        ~ G(0.001, 0.001)
  # tau.cst                      ~ G(0.001, 0.001)
  results.c_stwt_w <- run.jags("inst/individual_models/c_stwt_w.R",
                               data = data.jags,
                               n.chains = 2) # 832.0 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clinical.effect + wardtype.effect + sampletype.effect +
  #         ward.effect
  # clinical.effect   ~ N(0, 0.001)
  # wardtype.effect   ~ N(0, tau.wt)
  # sampletype.effect ~ N(0, tau.st)
  # ward.effect       ~ N(0, tau.w)
  # intercept         ~ N(0, 0.001)
  # tau.wt            ~ G(0.001, 0.001)
  # tau.st            ~ G(0.001, 0.001)
  # tau.w             ~ G(0.001, 0.001)
  results.cwtstw <- run.jags("inst/individual_models/cwtstw.R",
                             data = data.jags,
                             n.chains = 2) # 831.9 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clinical.effect + wardtype.effect + sampletype.effect +
  #         ward.effect + age.sq.effect * AGE^2_i + age.effect * AGE_i
  # clinical.effect   ~ N(0, 0.001)
  # wardtype.effect   ~ N(0, tau.wt)
  # sampletype.effect ~ N(0, tau.st)
  # ward.effect       ~ N(0, tau.w)
  # age.effect    ~ N(0, 0.001)
  # age.sq.effect ~ N(0, 0.001)
  # intercept         ~ N(0, 0.001)
  # tau.wt            ~ G(0.001, 0.001)
  # tau.st            ~ G(0.001, 0.001)
  # tau.w             ~ G(0.001, 0.001)
  results.cwtstwagesq <- run.jags("inst/individual_models/cwtstwagesq.R",
                                  data = data.jags,
                                  n.chains = 2) # 828.2 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_{ward,i})
  # p_{ward,i} = logit^{-1}(y_{ward,i})
  # y_{ward,i} = cstwt.effect_{ward}
  # cstwt.effect_{ward} ~ dnorm(cstwt.effect, tau.w)
  # cstwt.effect = intercept + clin.effect + wt.effect + st.effect
  # clin.effect   ~ dnorm(0, 0.001)
  # wt.effect     ~ N(0, tau.wt)
  # st.effect     ~ N(0, tau.st)
  # intercept     ~ N(0, 0.001)
  # tau.wt        ~ G(0.001, 0.001)
  # tau.st        ~ G(0.001, 0.001)
  # tau.w         ~ G(0.001, 0.001)
  results.cwtst_w <- run.jags("inst/individual_models/cwtst_w.R",
                              data = data.jags,
                              n.chains = 2) # 840.5 (sonia)

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = cstwt.effect_{ward} + age.sq.effect * age^2_i + age.effect * age_i
  # cstwt.effect_{ward} ~ dnorm(cstwt.effect, tau.w)
  # cstwt.effect = intercept + clin.effect + wt.effect + st.effect
  # clin.effect   ~ dnorm(0, 0.001)
  # wt.effect     ~ N(0, tau.wt)
  # st.effect     ~ N(0, tau.st)
  # age.effect    ~ N(0, 0.001)
  # age.sq.effect ~ N(0, 0.001)
  # intercept     ~ N(0, 0.001)
  # tau.wt        ~ G(0.001, 0.001)
  # tau.st        ~ G(0.001, 0.001)
  # tau.w         ~ G(0.001, 0.001)
  results.cwtst_wagesq <- run.jags("inst/individual_models/cwtst_wagesq.R",
                                   data = data.jags,
                                   n.chains = 2) # 837.8 (sonia) - takes ages

  # RESISTANCE_i ~ Bernoulli(p_i)
  # p_i = logit^{-1}(y_i)
  # y_i = intercept + clin.effect + wardtype.effect
  # clin.effect ~  dnorm(0, 0.001)
  # wardtype.effect ~  dnorm(0, tau.wt)
  # intercept     ~ N(0, 0.001)
  # tau.wt        ~ G(0.001, 0.001)
  results.cwt <- run.jags("inst/individual_models/cwt.R",
                          data = data.jags,
                          n.chains = 2) # 859.3 (sonia)

  # what are we looking for from the tau posterior? ----

  # sseff and psrf are calculated from monitored variables, so
  # then which variables do we want to monitor? ----

  models <- list(null = results.null,
                 g = results.g,
                 agegroup = results.agegroup,
                 agegroup2 = results.agegroup2,
                 age = results.age,
                 agesq = results.agesq,
                 agesqg = results.agesqg,
                 agesq_g = results.agesq_g,
                 c = results.c,
                 c_st = results.c_st,
                 cst = results.cst,
                 st = results.st,
                 h = results.h,
                 wt = results.wt,
                 h_w = results.h_w,
                 h_wt = results.h_wt,
                 h_st = results.h_st,
                 w = results.w,
                 wt_w = results.wt_w,
                 w_c = results.w_c,
                 w_st = results.w_st,
                 c_w = results.c_w,
                 wc = results.wc,
                 c_wt = results.c_wt,
                 w_wt_w = results.c_wt_w,
                 c_wt_st = results.c_wt_st,
                 c_wt_st_w = results.c_wt_st_w,
                 cwt_w = results.cwtw,
                 cwt_w = results.cwt_w,
                 cwt_w2 = results.cwt_w2,
                 cwtg_w = results.cwtg_w,
                 c_stwt_w = results.c_stwt_w,
                 cwtstw = results.cwtstw,
                 cwtstwagesq = results.cwtstwagesq,
                 cwtst_w = results.cwtst_w)

  tag <- gsub("/", "-", antibiotic_classes[x])
  saveRDS(models, paste0("inst/individual_results/", tag, "_results.rds"))

}


models <- dir("inst/individual_results", full.names = T) %>%
  .[grepl(".rds", .)] %>%
  lapply(function(x) readRDS(x))
names(models) <- gsub(".rds", "", dir("inst/individual_results")) %>%
  .[grep("_results",.)] %>%
  gsub("_results", "", .)


# Extract results from model objects --------------------------------------
results <- lapply(seq_along(models), function(x) {
  this_class <- models[[x]]
  model_names <- names(this_class)

  if(length(this_class) == 1 && is.na(this_class))
    return(data.frame(parameter = NA, Lower95 = NA, Median = NA,
                      Upper95 = NA, Mean = NA, SD = NA, Mode = NA,
                      MCerr = NA, MC.ofSD = NA, SSeff = NA, AC.10 = NA,
                      psrf = NA, pD = NA, DIC = NA, wAIC = NA, PSIS_loo = NA,
                      model = NA, antibiotic_class = antibiotic_classes[x]))

  tmp <- lapply(seq_along(this_class), function(y) {
    fits <- check_fit(this_class[[y]])

    inner_tmp <- this_class[[y]]$summaries %>%
      data.frame()

    inner_tmp <- cbind.data.frame(inner_tmp,
                                  parameter = rownames(inner_tmp)) %>%
      mutate(pD = mean(this_class[[y]]$pd),
             DIC = this_class[[y]]$dic$dic,
             wAIC = fits['waic'],
             PSIS_loo = fits['psis_loo'],
             model = model_names[y],
             antibiotic_class = antibiotic_classes[x])
  }) %>%
    do.call(rbind.data.frame, .)
}) %>%
  do.call(rbind.data.frame, .)




# Good psrf and effective sample size suggest good convergence, which can
# also be examined visually with a trace plot

# Check psrf, a measure of how well chains mix (the factor by which the
# variance in the estimate might be reduced with longer chains)
bad_psrf_values <- lapply(seq_along(models), function(x) {
  this.class <- models[[x]]
  lapply(seq_along(this.class), function(y){
    inner_tmp <- data.frame(this.class[[y]]$psrf$psrf) %>%
      select(.data$Point.est.)
    cbind.data.frame(inner_tmp, Parameter = rownames(inner_tmp)) %>%
      mutate(model = names(this.class)[y],
             antibiotic_class = names(models)[x])
    }) %>%
    do.call(rbind.data.frame, .) %>%
    # filter(Point.est. > 1.05)
    filter(Point.est. > 1.1)
}) %>%
  do.call(rbind.data.frame, .)
bad_psrf_values
write.csv(bad_psrf_values, "inst/individual_results/bad_psrf_values.csv")
# if DIC are comparable then find the best psrf






# Check effective sample size (equivalent number of independent interactions
# that the chain represents), which essentially accounts for autocorrelation
# in the chain (the correlation of estimates between one draw and the next).
# We want this to equal the number of posterior draws requested. In the
# presence of essentially no autocorrelation the values would be equal.
bad_sample_sizes <- lapply(seq_along(models), function(x) {
  this.class <- models[[x]]
  lapply(seq_along(this.class), function(y) {
    inner_tmp <- data.frame(SSeff = this.class[[y]]$mcse$sseff)
    cbind.data.frame(inner_tmp, Parameter = rownames(inner_tmp)) %>%
      mutate(model = names(this.class)[y],
             antibiotic_class = names(models)[x])
  }) %>%
    do.call(rbind.data.frame, .) %>%
    # filter(SSeff < 1000)
    filter(SSeff < 300)
}) %>%
  do.call(rbind.data.frame, .)
bad_sample_sizes
write.csv(bad_sample_sizes, "individual_results/bad_sample_sizes.csv")
#extend.jags





# Sort models by DIC for each antibiotic class ----------------------------

# DIC can be unreliable for complex hierarchical models, so use caution when
# applying DIC to model selection problems. The WAIC metric is generally
# considered a better metric than DIC- it is applicable more widely than
# DIC and has a more solid theoretical foundation.
best_models <- results %>%
  select(antibiotic_class, model, DIC) %>%
  unique() %>%
  arrange(antibiotic_class, DIC) %>%
  group_by(antibiotic_class) %>%
  mutate(index = row_number()) %>%
  reshape2::dcast(index ~ antibiotic_class, value.var = "model") %>%
  select(-index)
write.csv(best_models, "individual_results/models_DIC.csv")

plot_this <- results %>%
  select(antibiotic_class, model, DIC) %>%
  unique() %>%
  arrange(antibiotic_class, DIC)  %>%
  group_by(antibiotic_class) %>%
  mutate(index = row_number(),
         DIC = DIC - min(DIC))

best_ones <- plot_this %>%
  filter(DIC < 6) %>%
  group_by(model) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
write.csv(best_ones, "individual_results/best_ones.csv")


p <- merge(plot_this, bad_sample_sizes, all.x = T) %>%
  merge(bad_psrf_values, all.x = T) %>%
  rename(psrf = Point.est.) %>%
  mutate(fill = case_when(DIC < 2 ~ "<2",
                          DIC < 4 ~ "<4",
                          DIC < 6 ~ "<6",
                          DIC < 8 ~ "<8",
                          DIC < 10 ~ "<10",
                          T ~ "bad"),
         psrf_tag = case_when(is.na(psrf) ~ "good",
                              T ~ "bad"),
         sseff_tag = case_when(is.na(SSeff) ~ "good",
                               T ~ "bad")) %>%
  mutate(fill = factor(fill, levels = c("<2", "<4", "<6", "<8", "<10", "bad"))) %>%
  group_by(antibiotic_class, model) %>%
  mutate(bad = case_when(
    any(psrf_tag == "bad") & any(sseff_tag == "bad") ~ paste0(model, "^p*{}['s']"),
    any(psrf_tag == "bad") ~ paste0(model, "^p"),
    any(sseff_tag == "bad") ~ paste0(model, "[s]"),
    T ~ model)) %>%
  select(antibiotic_class, model, index, fill, bad) %>%
  unique()


g <- ggplot2::ggplot(p) + theme_minimal() +
  # ggplot2::scale_fill_gradientn(
  #   colours = c('#ffffcc', '#fed976', 'white',
  #               '#74c476','#00441b', 'white',
  #               '#e7298a','#67001f'),
  #   values = c(0, 1.9, 2,
  #              2.1, 9.9, 10,
  #              10.1, max(plot_this$DIC, na.rm = T))/100) +
  ggplot2::geom_tile(ggplot2::aes(x = index, y = antibiotic_class,
                                  group = antibiotic_class,
                                  fill = fill)) +
  ggplot2::scale_fill_manual(values = c("#ffff99", "#edf8fb", "#b2e2e2",
                                        "#66c2a4", "#238b45", "#bababa")) +
  ggplot2::geom_text(ggplot2::aes(x = index, y = antibiotic_class,
                                  label = bad),
                     parse = T) +
  ggplot2::labs(y = "Antibiotic class", fill = "DIC") +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())

ggsave("individual_results/models_DIC.pdf", g, width = 40, height = 8)





# Sort models by wAIC for each antibiotic class ---------------------------

# best_models <- results %>%
#   select(antibiotic_class, model, wAIC) %>%
#   unique() %>%
#   arrange(antibiotic_class, wAIC) %>%
#   group_by(antibiotic_class) %>%
#   mutate(index = row_number()) %>%
#   reshape2::dcast(index ~ antibiotic_class, value.var = "model") %>%
#   select(-index)
# write.csv(best_models, "models_wAIC.csv")
#
# plot_this <- results %>%
#   select(antibiotic_class, model, wAIC) %>%
#   unique() %>%
#   arrange(antibiotic_class, wAIC) %>%
#   group_by(antibiotic_class) %>%
#   mutate(index = row_number(),
#          wAIC = wAIC - min(wAIC))
#
# g <- plot_this %>% ggplot2::ggplot() + theme_minimal() +
#   ggplot2::scale_fill_gradientn(
#     colours = c('#ffffcc', '#fed976', 'white',
#                 '#74c476','#00441b', 'white',
#                 '#e7298a','#67001f'),
#     values = c(0, 1.9, 2,
#                2.1, 9.9, 10,
#                10.1, max(plot_this$wAIC, na.rm = T))/100) +
#   ggplot2::geom_tile(ggplot2::aes(x = index, y = antibiotic_class,
#                                   group = antibiotic_class,
#                                   fill = wAIC)) +
#   ggplot2::geom_text(ggplot2::aes(x = index, y = antibiotic_class,
#                                   label = model))
#
# ggsave("models_wAIC.pdf", g, width = 20, height = 8)




# Sort models by PSIS-loo for each antibiotic class -----------------------

# best_models <- results %>%
#   select(antibiotic_class, model, PSIS_loo) %>%
#   unique() %>%
#   arrange(antibiotic_class, PSIS_loo) %>%
#   group_by(antibiotic_class) %>%
#   mutate(index = row_number()) %>%
#   reshape2::dcast(index ~ antibiotic_class, value.var = "model") %>%
#   select(-index)
# write.csv(best_models, "models_PSIS_loo.csv")
#
# plot_this <- results %>%
#   select(antibiotic_class, model, PSIS_loo) %>%
#   unique() %>%
#   arrange(antibiotic_class, PSIS_loo) %>%
#   group_by(antibiotic_class) %>%
#   mutate(index = row_number(),
#          PSIS_loo = PSIS_loo - min(PSIS_loo))
#
# g <- plot_this %>% ggplot2::ggplot() + theme_minimal() +
#   ggplot2::scale_fill_gradientn(
#     colours = c('#ffffcc', '#fed976', 'white',
#                 '#74c476','#00441b', 'white',
#                 '#e7298a','#67001f'),
#     values = c(0, 1.9, 2,
#                2.1, 9.9, 10,
#                10.1, max(plot_this$PSIS_loo, na.rm = T))/100) +
#   ggplot2::geom_tile(ggplot2::aes(x = index, y = antibiotic_class,
#                                   group = antibiotic_class,
#                                   fill = PSIS_loo)) +
#   ggplot2::geom_text(ggplot2::aes(x = index, y = antibiotic_class,
#                                   label = model))
#
# ggsave("models_PSIS_loo.pdf", g, width = 20, height = 8)










# Monte Carlo error is an estimate of the uncertainty contributed by
# only having a finite number of posterior draws. We want less than 5% of the
# posterior standard deviation.
# model$mcse
bad_mcerr <- lapply(seq_along(plot_these), function(x) {
  this.class <- plot_these[[x]]
  lapply(seq_along(this.class), function(y) {
    inner_tmp <- data.frame(mcerr = this.class[[y]]$mcse$mcse /
                 this.class[[y]]$mcse$ssd * 100)
     cbind.data.frame(inner_tmp, Parameter = rownames(inner_tmp)) %>%
      mutate(model = names(this.class)[y],
             antibiotic_class = names(plot_these)[x])
  }) %>%
    do.call(rbind.data.frame, .) %>%
    filter(mcerr < 5)
}) %>%
  do.call(rbind.data.frame, .)
bad_mcerr
write.csv(bad_mcerr, "individual_results/bad_mcerr.csv")









# # Generate lookup table for plots
# lookup <- get_data(classification = antibiotic_classes[2],
#                    kleb = "Klebsiella pneumoniae")$lookup_tables
#
# # Plot clinical
# g <- lapply(seq_along(plot_these), function(x) {
#   label <- data.frame(Parameter = paste0("c.effect[",
#                                          lookup$clinical$index, "]"),
#                       Label = lookup$clinical$clinical) %>%
#     mutate(Label = if_else(Label == "yes", "clinical", "carriage"))
#   model <- plot_these[[x]]$clinical
#   p <- get_parameters(model) %>%
#     .[grepl("effect", .)] %>%
#     plot_jags(model, ., label, 0)
#
#   title <- cowplot::ggdraw() + cowplot::draw_label(names(plot_these)[x],
#                                                    fontface = "bold")
#   cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
# }) %>%
#   cowplot::plot_grid(plotlist = ., ncol = 2)
# ggsave("clinical.pdf", g, height = 16, width = 10)
#
#
#
#
# # Extract the best model for each antibiotic class and convert log-odds
# # predictions to odds and probabilities
# odds <- lapply(seq_along(plot_these), function(x) {
#   tmp <- best_models %>% select(names(plot_these)[x]) %>% .[1,1]
#   this.model <- plot_these[[x]][[tmp]]
#
#   data.frame(this.model$summaries)[, "Mean", drop = F] %>%
#     tibble::rownames_to_column("variable") %>%
#     mutate(odds = case_when(
#       grepl("effect", variable) ~ exp(Mean),
#       grepl("diff", variable) ~ exp(Mean)),
#       probability = case_when(
#         grepl("effect", variable) ~ plogis(Mean),
#         grepl("diff", variable) ~ plogis(Mean)))
# })
# names(odds) <- names(plot_these)
# odds
#
#
#
#
#
#
# # dnorm(0, 0.001)
# qplot(rnorm(1e5, mean = 0, sd = sqrt(1/0.001)))



# Check autocorrelation of best models?
# Should I use autorun.jags?






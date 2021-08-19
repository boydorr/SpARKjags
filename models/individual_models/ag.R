model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability agegroup.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(agegroup.prob[age_group[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(agegroup.prob[age_group[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(agegroup.prob[age_group[v_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(agegroup.prob[age_group[o_sample_GUID[o]]])
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for agegroup.effect (log-odds for each age group). Sample
  # different agegroup.effect from normal distribution for each age group and
  # convert to a probability). Since there is only one response variable, put
  # intercept here.
  for (a in 1:N_age_group)
  {
    agegroup.effect[a] ~ dnorm(intercept, tau.agegroup) # like rnorm in R
    logit(agegroup.prob[a]) <- agegroup.effect[a]
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.0001)

  # Prior values for precision
  tau.agegroup ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.agegroup <- sqrt(1/tau.agegroup)

  #monitor# full.pd, dic, deviance, agegroup.effect, agegroup.prob, intercept, sd.agegroup
}

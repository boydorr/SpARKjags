model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability age.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(age.prob[h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(age.prob[gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(age.prob[v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(age.prob[o_sample_GUID[o]])
  }

  # ------------------------

  # Prior distribution for age.effect (since there is only one response
  # variable, put intercept here; putting the intercept here rather than
  # above stops variables from covarying negatively)
  for (s in 1:N_sample)
  {
    logit(age.prob[s]) <- intercept + age.effect * age[s]
  }

  # Prior distributions for estimates (because age is continuous and has a
  # single slope, precision takes a single value)
  age.effect ~ dnorm(0, 0.001)

  # ------------------------

  intercept ~ dnorm(0, 0.001)

  #monitor# full.pd, dic, deviance, age.effect, intercept
}

model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wt.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wt.prob[ward_type[ward[h_sample_GUID[p]]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gp.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(v.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(o.prob)
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Since there is only one response variable, put
  # intercept here.
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  gp.effect ~ dnorm(intercept, tau.wt) # equivalent to wt.effect
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(intercept, tau.wt)  # equivalent to wt.effect
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(intercept, tau.wt)  # equivalent to wt.effect
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- sqrt(1/tau.wt)

  #monitor# full.pd, dic, deviance, gp.prob, v.prob, o.prob, intercept, sd.wt
}

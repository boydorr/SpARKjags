model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability w.prob (gp.prob, v.prob,
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

  # Prior distribution for h.effect (log-odds for each hospital). Sample
  # different h.effect from normal distribution for each hospital and
  # convert to a probability).
  for (h in 1:N_hosp)
  {
    h.effect[h] ~ dnorm(intercept, tau.hosp)
    logit(h.prob[h]) <- h.effect[h]
  }

  # equivalent to h.effect
  nh.effect ~ dnorm(intercept, tau.hosp)
  logit(nh.prob) <- nh.effect

  # ------------------------

  # Define the priors:
  # Prior distribution for ward.effect (log-odds for each hospital ward). Sample
  # different ward.effect from normal distribution for each hospital ward and
  # convert to a probability).
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(h.effect[hospital[wt]], tau.wt)
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  # equivalent to ward.effect
  gp.effect ~ dnorm(nh.effect, tau.wt)
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(nh.effect, tau.wt)
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(nh.effect, tau.wt)
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.hosp ~ dgamma(0.001, 0.001)
  tau.wt ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.hosp <- sqrt(1/tau.hosp)
  sd.wt <- sqrt(1/tau.wt)

  #monitor# full.pd, dic, deviance, intercept, h.prob, nh.prob, sd.hosp, sd.wt
}

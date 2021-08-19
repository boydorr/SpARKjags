model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability ward.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(ward.prob[ward[h_sample_GUID[p]]])
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

  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Assuming the presence of resistance in a
  # hospital ward is based on ward type there, put intercept here.
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
  }

  # equivalent to wt.effect
  nh.effect ~ dnorm(intercept, tau.wt)
  logit(nh.prob) <- nh.effect



  # Prior distribution for wtw.effect (log-odds for each hospital ward). Sample
  # different wtw.effect from normal distribution for each hospital ward and
  # convert to a probability).
  for (w in hosp_wards)
  {
    wtw.effect[w] ~ dnorm(wt.effect[ward_type[w]], tau.ward)
    logit(ward.prob[w]) <- wtw.effect[w]
  }

  # equivalent to wtw.effect
  gp.effect ~ dnorm(nh.effect, tau.ward)
  v.effect ~ dnorm(nh.effect, tau.ward)
  o.effect ~ dnorm(nh.effect, tau.ward)

  # convert to probability
  logit(gp.prob) <- gp.effect
  logit(v.prob) <- v.effect
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.ward ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- sqrt(1/tau.wt)
  sd.ward <- sqrt(1/tau.ward)

  #monitor# full.pd, dic, deviance, intercept, nh.prob, sd.wt, sd.ward
}

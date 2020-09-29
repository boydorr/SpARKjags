model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wc.prob[ward[h_sample_GUID[p]],
                                clinical[sample_type[h_sample_GUID[p]]]])
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
  # Prior distribution for ward.effect (log-odds for each hospital ward). Sample
  # different ward.effect from normal distribution for each hospital ward and
  # convert to a probability). Since there is only one response variable, put
  # intercept here.
  for (w in hosp_wards)
  {
    ward.effect[w] ~ dnorm(intercept, tau.ward)

    for (c in c(ncarr, nclin))
    {
      wc.effect[w,c] ~ dnorm(ward.effect[w], tau.wc)
      logit(wc.prob[w,c]) <- wc.effect[w,c]
    }
  }

  # equivalent to ward.effect
  nh.effect ~ dnorm(intercept, tau.ward)
  logit(nh.prob) <- nh.effect

  # equivalent to wc.effect
  gp.effect ~ dnorm(nh.effect, tau.wc)
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(nh.effect, tau.wc)
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(nh.effect, tau.wc)
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.ward ~ dgamma(0.001, 0.001)
  tau.wc ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.ward <- 1/sqrt(tau.ward)
  sd.wc <- 1/sqrt(tau.wc)

  #monitor# full.pd, dic, deviance, gp.prob, v.prob, o.prob, sd.ward, sd.sc
}

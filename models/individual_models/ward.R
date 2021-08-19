model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability ward.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(ward.prob[ward[h_sample_GUID[p]]])
    # For WAIC computation
    # h_like[p] <- logdensity.bin(h_resist[p],
    #                             ward.prob[ward[h_sample_GUID[p]]], 1)
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gp.prob)
    # For WAIC computation
    # gp_like[gp] <- logdensity.bin(gp_resist[gp], gp.prob, 1)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(v.prob)
    # For WAIC computation
    # v_like[v] <- logdensity.bin(v_resist[v], v.prob, 1)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(o.prob)
    # For WAIC computation
    # o_like[o] <- logdensity.bin(o_resist[o], o.prob, 1)
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for ward.effect (log-odds for each hospital ward). Sample
  # different ward.effect from normal distribution for each hospital ward and
  # convert to a probability). Since there is only one response variable, put
  # intercept here.
  for (w in hosp_wards)
  {
    ward.effect[w] ~ dnorm(intercept, tau.ward) # params = 119 + 1 + 1
    logit(ward.prob[w]) <- ward.effect[w]
  }

  # equivalent to ward.effect
  gp.effect ~ dnorm(intercept, tau.ward)        # params = 1 + 1 + 1
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(intercept, tau.ward)         # params = 1 + 1 + 1
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(intercept, tau.ward)         # params = 1 + 1 + 1
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)                   # params = 1

  # Prior values for precision
  tau.ward ~ dgamma(0.001, 0.001)               # params = 1

  # Convert precisions to sd
  sd.ward <- sqrt(1/tau.ward)

  #monitor# full.pd, dic, deviance, intercept, gp.prob, v.prob, o.prob, sd.ward
}

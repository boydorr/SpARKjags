model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability h.prob (or nh.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(h.prob[hospital[ward[h_sample_GUID[p]]]])
    # For WAIC computation
    # h_like[p] <- logdensity.bin(h_resist[p],
    #                             h.prob[hospital[ward[h_sample_GUID[p]]]], 1)
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(nh.prob)
    # For WAIC computation
    # gp_like[gp] <- logdensity.bin(gp_resist[gp], nh.prob, 1)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(nh.prob)
    # For WAIC computation
    # v_like[v] <- logdensity.bin(v_resist[v], nh.prob, 1)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(nh.prob)
    # For WAIC computation
    # o_like[o] <- logdensity.bin(o_resist[o], nh.prob, 1)
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for h.effect (log-odds for each hospital). Sample
  # different h.effect from normal distribution for each hospital and
  # convert to a probability). Since there is only one explanatory variable,
  # put intercept here.
  for (h in 1:N_hosp)
  {
    h.effect[h] ~ dnorm(intercept, tau)
    logit(h.prob[h]) <- h.effect[h]
  }

  # nh is considered a hospital
  nh.effect ~ dnorm(intercept, tau)
  logit(nh.prob) <- nh.effect

  # ------------------------

  # Prior value for intercept (log-odds of the average resistance in all samples)
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)

  #monitor# full.pd, dic, deviance, intercept, h.prob, nh.prob, sd
}

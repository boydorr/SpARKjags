model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability h.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(h.prob)
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
  logit(h.prob) <- intercept
  logit(o.prob) <- intercept
  logit(gp.prob) <- intercept
  logit(v.prob) <- intercept

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, tau)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)

  #monitor# full.pd, dic, deviance, h.prob, intercept, sd
}

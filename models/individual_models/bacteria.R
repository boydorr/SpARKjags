model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability b.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(b_prob[h_bacteria[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(b_prob[gp_bacteria[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(b_prob[v_bacteria[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(b_prob[o_bacteria[o]])
  }

  # ------------------------

  for (b in bact_species)
  {
    logit(b_prob[b]) <- b_value[b]
    b_value[b] ~ dnorm(mu, tau)
  }

  # ------------------------

  mu ~ dnorm(0, 0.0001)
  tau ~ dgamma(0.001, 0.001)

  #monitor# full.pd, dic, deviance, tau, mu, b_prob
}

model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability g.prob
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(g.prob[gender[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(g.prob[gender[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(g.prob[gender[v_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(g.prob[gender[o_sample_GUID[o]]])
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for gender.effect (log-odds for each gender). Sample
  # different gender.effect from normal distribution for each gender and
  # convert to a probability). Since there is only one explanatory variable,
  # with only two elements, set intercept as 0, because if you have more
  # degrees of freedom than levels in the variable, then the caterpillar might
  # have drifting issues.
  for (g in genders)
  {
    gender.effect[g] ~ dnorm(0, tau)
    logit(g.prob[g]) <- gender.effect[g]
  }

  # ------------------------

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- sqrt(1/tau)

  # Calculate odds
  g.diff <- gender.effect[male] - gender.effect[female]
  odds.g <- exp(g.diff)

  #monitor# full.pd, dic, deviance, gender.effect, g.prob, g.diff, odds.g, sd
}

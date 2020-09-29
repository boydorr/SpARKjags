model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wtc.prob[ward_type[ward[h_sample_GUID[p]]],
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
  # Prior distribution for clin.effect (log-odds for each clinical class).
  # Sample different clin.effect from normal distribution for each clinical
  # class and convert to a probability). Since there is only one response
  # variable, with only two levels, set intercept to 0.
  #
  # Prior distribution for cwt.effect (log-odds for each ward type in each
  # clinical class). Sample different cw.effect from normal distribution (with
  # mean clin.effect) for each ward type in each clinical class and convert
  # to a probability.

  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]


  # two level hierarchical model

  for (c in c(ncarr, nclin)) # 1, 2!
  {
    for (wt in hosp_wardtypes)
    {
      # clinical effect is different for each ward type
      cwt.effect[wt,c] ~ dnorm(clin.effect[c], tau)
      logit(wtc.prob[wt,c]) <- cwt.effect[wt,c]
    }
  }

  # equivalent to cwt.effect
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau)
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(clin.effect[v_clinical], tau)
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(clin.effect[o_clinical], tau)
  logit(o.prob) <- o.effect

  # ------------------------

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- 1/sqrt(tau)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, clin.effect, gp.prob, v.prob, o.prob, odds.c, sd
}

model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(stc.prob[sample_type[h_sample_GUID[p]],
                                 clinical[sample_type[h_sample_GUID[p]]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(stc.prob[sample_type[gp_sample_GUID[gp]],
                                   gp_clinical])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(stc.prob[sample_type[v_sample_GUID[v]],
                                 v_clinical])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(stc.prob[sample_type[o_sample_GUID[o]],
                                 o_clinical])
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Put intercept here.

  # Prior distribution for wc.effect (log-odds for each clinical class). Sample
  # different clin.effect from normal distribution for each clinical class and
  # convert to a probability). Since intercept is in wt, set mean to 0.

  # Prior distribution for wtc.effect
  for (c in c(ncarr, nclin)) # 1, 2!
  {
    clin.effect[c] ~ dnorm(0, tau.clin)
    logit(c.prob[c]) <- clin.effect[c]

    for (st in sampletypes)
    {
      stc.effect[st,c] ~ dnorm(clin.effect[c], tau.stc)
      logit(stc.prob[st,c]) <- stc.effect[st,c]
    }
  }

  # ------------------------

  # Prior values for precision
  tau.clin ~ dgamma(0.001, 0.001)
  tau.stc ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.clin <- 1/sqrt(tau.clin)
  sd.stc <- 1/sqrt(tau.stc)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, clin.effect, c.prob, odds.c, sd.clin, sd.stc
}

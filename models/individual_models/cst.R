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
  # Prior distribution for stc.effect (log-odds for each clinical class). Sample
  # different clin.effect from normal distribution for each clinical class and
  # convert to a probability). Since intercept is in stc, set mean to 0.
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  # Prior distribution for w.effect
  for (st in sampletypes)
  {
    st.effect[st] ~ dnorm(0, tau)
  }

  # Prior distribution for stc.effect
  for (c in c(ncarr, nclin)) # 1, 2!
  {
    for (st in sampletypes)
    {
      stc.effect[st,c] <- intercept + st.effect[st] + clin.effect[c]
      logit(stc.prob[st,c]) <- stc.effect[st,c]
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd <- 1/sqrt(tau)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, clin.effect, odds.c, c.diff, sd
}

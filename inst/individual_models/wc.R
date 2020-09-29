model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability cw.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cw.prob[ward[h_sample_GUID[p]],
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
  # class and convert to a probability).
  #
  # Since we're estimating exactly one parameter (the carriage effect),
  # we can't use more than one thing to explain it... so we don't include a
  # hyperparameter (e.g. mu or tau).
  #
  # log odds of a carriage sample = intercept + clin.effect[ncarr]
  # log odds of a clinical sample = intercept - clin.effect[ncarr]
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  # Prior distribution for ward.effect
  for (w in hosp_wards)
  {
    # for a random effect, report the estimated variance between wards
    ward.effect[w] ~ dnorm(0, tau)
  }

  for (c in c(ncarr, nclin))
  {
    for (w in hosp_wards)
    {
      cw.effect[w,c] <- intercept + ward.effect[w] + clin.effect[c]
      logit(cw.prob[w,c]) <- cw.effect[w,c]
    }
  }

  # equivalent to clin.effect + ward.effect? intercept? tau? ----!!!!
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau)
  logit(gp.prob) <- gp.effect

  v.effect ~ dnorm(clin.effect[v_clinical], tau)
  logit(v.prob) <- v.effect

  o.effect ~ dnorm(clin.effect[o_clinical], tau)
  logit(o.prob) <- o.effect

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

  #monitor# full.pd, dic, deviance, intercept, clin.effect, gp.prob, v.prob, o.prob, odds.c, c.diff, sd
}

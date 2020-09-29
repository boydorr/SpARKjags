model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability cw.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cwtw.prob[clinical[sample_type[h_sample_GUID[p]]],
                                  ward_type[ward[h_sample_GUID[p]]],
                                  ward[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gpw.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(vw.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(ow.prob)
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
    ward.effect[w] ~ dnorm(0, tau.w)
  }

  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
  }

  for (c in c(ncarr, nclin))
  {
    for (wt in hosp_wardtypes)
    {
      for (w in hosp_wards)
      {
        cwtw.effect[c,wt,w] <- intercept + clin.effect[c] + wt.effect[wt] +
          ward.effect[w]
        logit(cwtw.prob[c,wt,w]) <- cwtw.effect[c,wt,w]
      }
    }
  }

  # equivalent to clin + wt
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau.wt)
  v.effect ~ dnorm(clin.effect[v_clinical], tau.wt)
  o.effect ~ dnorm(clin.effect[o_clinical], tau.wt)

  # equivalent to clin + wt + w
  gpw.effect ~ dnorm(gp.effect, tau.w)
  vw.effect ~ dnorm(v.effect, tau.w)
  ow.effect ~ dnorm(o.effect, tau.w)


  # convert to probability
  logit(gpw.prob) <- gpw.effect
  logit(vw.prob) <- vw.effect
  logit(ow.prob) <- ow.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.w ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.w <- 1/sqrt(tau.w)
  sd.wt <- 1/sqrt(tau.wt)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, clin.effect, gp.prob, v.prob, o.prob, odds.c, c.diff, sd.w, sd.wt
}

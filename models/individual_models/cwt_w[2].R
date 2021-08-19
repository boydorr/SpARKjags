model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wc.prob[ward[h_sample_GUID[p]],
                                clinical[sample_type[h_sample_GUID[p]]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(wgp.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(wv.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(wo.prob)
  }

  # ------------------------

  # Define the priors:
  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Put intercept here.
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  # Prior distribution for wc.effect (log-odds for each clinical class). Sample
  # different clin.effect from normal distribution for each clinical class and
  # convert to a probability). Since intercept is in wt, set mean to 0.
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  # Prior distribution for wtc.effect
  for (c in c(ncarr, nclin)) # 1, 2!
  {
    # Ward type and clinical state are independent of each other
    for (wt in hosp_wardtypes)
    {
      wtc.effect[wt,c] <- clin.effect[c] + wt.effect[wt]
    }

    # Hospital ward (wc) is nested within ward type & clinical state
    for (w in hosp_wards)
    {
      wc.effect[w,c] ~ dnorm(wtc.effect[ward_type[w],c], tau.w)
      logit(wc.prob[w,c]) <- wc.effect[w,c]
    }

  }

  # equivalent to wtc.effect
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau.wt)
  v.effect ~ dnorm(clin.effect[v_clinical], tau.wt)
  o.effect ~ dnorm(clin.effect[o_clinical], tau.wt)

  # equivalent to wc.effect
  wgp.effect ~ dnorm(gp.effect, tau.w)
  wv.effect ~ dnorm(v.effect, tau.w)
  wo.effect ~ dnorm(o.effect, tau.w)

  # convert to probability
  logit(wgp.prob) <- wgp.effect
  logit(wv.prob) <- wv.effect
  logit(wo.prob) <- wo.effect

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.w ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- 1/sqrt(tau.wt)
  sd.w <- 1/sqrt(tau.w)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, clin.effect, wgp.prob, wv.prob, wo.prob, odds.c, c.diff, sd.wt, sd.w
}

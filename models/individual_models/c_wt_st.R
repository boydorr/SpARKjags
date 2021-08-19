model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cwtst.prob[clinical[sample_type[h_sample_GUID[p]]],
                                  ward_type[ward[h_sample_GUID[p]]],
                                  sample_type[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gpst.prob)
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(vst.prob)
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(ost.prob)
  }

  # ------------------------

  # Define the priors:
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  for (c in c(ncarr, nclin)) # 1, 2!
  {
    for (wt in hosp_wardtypes)
    {
      cwt.effect[c,wt] ~ dnorm(clin.effect[c], tau.wt)

      for (st in sampletypes)
      {
        cwtst.effect[c,wt,st] ~ dnorm(cwt.effect[c,wt], tau.st)
        logit(cwtst.prob[c,wt,st]) <- cwtst.effect[c,wt,st]
      }
    }
  }

  # equivalent to cwt.effect
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau.wt)
  v.effect ~ dnorm(clin.effect[v_clinical], tau.wt)
  o.effect ~ dnorm(clin.effect[o_clinical], tau.wt)

  # equivalent to cwtst.effect
  gpst.effect ~ dnorm(gp.effect, tau.st)
  vst.effect ~ dnorm(v.effect, tau.st)
  ost.effect ~ dnorm(o.effect, tau.st)

  # convert to probability
  logit(gpst.prob) <- gpst.effect
  logit(vst.prob) <- vst.effect
  logit(ost.prob) <- ost.effect

  # ------------------------

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.st ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- 1/sqrt(tau.wt)
  sd.st <- 1/sqrt(tau.st)

  # Calculate odds
  c.diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds.c <- exp(c.diff)

  #monitor# full.pd, dic, deviance, intercept, clin.effect, gp.prob, v.prob, o.prob, odds.c, c.diff, sd.wt, sd.st
}

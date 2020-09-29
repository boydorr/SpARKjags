model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cstwtw.prob[ward[h_sample_GUID[p]],
                                    clinical[sample_type[h_sample_GUID[p]]],
                                    sample_type[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gpcst.prob[gp_clinical,
                                     sample_type[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(vcst.prob[v_clinical,
                                  sample_type[v_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(ocst.prob[o_clinical,
                                  sample_type[o_sample_GUID[o]]])
  }

  # ------------------------

  # Define the priors:
  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]

  for (st in sampletypes)
  {
    st.effect[st] ~ dnorm(0, tau.st)
    logit(st.prob[st]) <- st.effect[st]
  }

  # ------------------------

  # Prior distribution for wt.effect (log-odds for each ward type). Sample
  # different wt.effect from normal distribution for each ward type and
  # convert to a probability). Put intercept here.
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  # equivalent to wt.effect
  nh.effect ~ dnorm(intercept, tau.wt)
  logit(nh.prob) <- nh.effect

  # ------------------------

  # Prior distribution for w.effect (log-odds for each ward). Sample
  # different ward.effect from normal distribution for each clinical class and
  # convert to a probability). w.effect depends on wardtype effect.
  for (w in hosp_wards)
  {
    w.effect[w] ~ dnorm(0, tau.w)
  }

  # equivalent to wt.effect + w.effect
  gp.effect ~ dnorm(nh.effect, tau.w)
  v.effect ~ dnorm(nh.effect, tau.w)
  o.effect ~ dnorm(nh.effect, tau.w)

  # ------------------------

  # Prior distribution for st.effect (log-odds for each sample type). Sample
  # different st.effect from normal distribution for each sample type and
  # convert to a probability). cst.effect depends on clin.effect and is
  # different for each sampletype. cst.effect (clinical state) and w.effect are
  # independent. cstwtw.effect is different for each ward, clinical state, and
  # sampletype.
  for (c in c(ncarr, nclin))
  {
    for (st in sampletypes)
    {
      cst.effect[c,st] <- clin.effect[c] + st.effect[st]

      for (w in hosp_wards)
      {
        cstwtw.effect[w,c,st] <- cst.effect[c,st] + wt.effect[ward_type[w]] +
          w.effect[w]
        logit(cstwtw.prob[w,c,st]) <- cstwtw.effect[w,c,st]
      }

      # equivalent to cstwtw.prob
      logit(gpcst.prob[c,st]) <- cst.effect[c,st] + gp.effect
      logit(vcst.prob[c,st]) <- cst.effect[c,st] + v.effect
      logit(ocst.prob[c,st]) <- cst.effect[c,st] + o.effect
    }
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)
  tau.w ~ dgamma(0.001, 0.001)
  tau.st ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.wt <- 1/sqrt(tau.wt)
  sd.w <- 1/sqrt(tau.w)
  sd.st <- 1/sqrt(tau.st)

  #monitor# full.pd, dic, deviance, intercept, gp.prob, v.prob, o.prob, sd.wt, sd.w
}

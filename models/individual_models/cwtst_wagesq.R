model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wc.prob (gp.prob, v.prob,
  # and o.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(cstwtw.prob[clinical[sample_type[h_sample_GUID[p]]],
                                    sample_type[h_sample_GUID[p]],
                                    ward[h_sample_GUID[p]],
                                    h_sample_GUID[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gp.prob[gp_clinical,
                                  sample_type[gp_sample_GUID[gp]],
                                  gp_sample_GUID[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(v.prob[v_clinical,
                               sample_type[v_sample_GUID[v]],
                               v_sample_GUID[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(o.prob[o_clinical,
                               sample_type[o_sample_GUID[o]],
                               o_sample_GUID[o]])
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

  # ------------------------

  for (s in 1:N_sample)
  {
    mu.age[s] <- age.effect * age[s] + age.sq.effect * pow(age[s], 2)
  }

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
      cstage2.effect[c,st] <- clin.effect[c] + st.effect[st]

      for (wt in hosp_wardtypes)
      {
        cstwt.effect[c,st,wt] <- cstage2.effect[c,st] + wt.effect[wt]
      }

      for (w in hosp_wards)
      {
        cstwtw.effect[c,st,w] ~ dnorm(cstwt.effect[c,st,ward_type[w]], tau.w)
      }

      # equivalent to cstwt.effect
      cstage2nh.effect[c,st] <- cstage2.effect[c,st] + nh.effect

      # equivalent to cstwtw.effect
      gp.effect[c,st] ~ dnorm(cstage2nh.effect[c,st], tau.w)
      v.effect[c,st] ~ dnorm(cstage2nh.effect[c,st], tau.w)
      o.effect[c,st] ~ dnorm(cstage2nh.effect[c,st], tau.w)
    }

    for (s in 1:N_sample)
    {
      for (w in hosp_wards)
      {
        logit(cstwtw.prob[c,sample_type[s],w,s]) <-
          cstwtw.effect[c,sample_type[s],w] + mu.age[s]
      }

      # convert to probability
      logit(gp.prob[c,sample_type[s],s]) <- gp.effect[c,sample_type[s]] + mu.age[s]
      logit(v.prob[c,sample_type[s],s]) <- v.effect[c,sample_type[s]] + mu.age[s]
      logit(o.prob[c,sample_type[s],s]) <- o.effect[c,sample_type[s]] + mu.age[s]

    }
  }
  # ------------------------

  age.effect ~ dnorm(0, 0.001)
  age.sq.effect ~ dnorm(0, 0.001)

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

  # Calculate odds
  c_diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds_c <- exp(c_diff)

  #monitor# full.pd, dic, deviance, c_diff, odds_c, sd.wt, sd.w, sd.st
}

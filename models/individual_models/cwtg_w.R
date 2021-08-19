model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wcg.prob (gpg.prob, vg.prob,
  # and og.prob)
  for (p in 1:N_patients)    # hospital samples
  {
    h_resist[p] ~ dbern(cwtg_w.prob[ward[h_sample_GUID[p]],
                                   clinical[sample_type[h_sample_GUID[p]]],
                                   gender[h_sample_GUID[p]]])
  }

  for (gp in 1:N_gp)         # gp samples
  {
    gp_resist[gp] ~ dbern(gpg.prob[gender[gp_sample_GUID[gp]]])
  }

  for (v in 1:N_volunteers)  # volunteer samples
  {
    v_resist[v] ~ dbern(vg.prob[gender[o_sample_GUID[v]]])
  }

  for (o in 1:N_outpatients) # outpatient samples
  {
    o_resist[o] ~ dbern(og.prob[gender[o_sample_GUID[o]]])
  }

  # ------------------------

  # Prior distribution for gender.effect (log-odds for each gender).
  # Sample different gender.effect from normal distribution for each gender
  # and convert to a probability). Since there is only one response
  # variable, with only two elements, set intercept as 0.
  gender.effect[female] ~ dnorm(0, 0.001)             # params = 1
  gender.effect[male] <- -gender.effect[female]

  # Prior distribution for wt.effect (log-odds for each ward type).
  # Sample different wt.effect from normal distribution for each ward type
  # and convert to a probability).
  for (wt in hosp_wardtypes)
  {
    wt.effect[wt] ~ dnorm(intercept, tau.wt)          # params = 8 + 1
    logit(wt.prob[wt]) <- wt.effect[wt]
  }

  # Prior distribution for clin.effect (log-odds for each clinical class).
  # Sample different clin.effect from normal distribution for each clinical
  # class and convert to a probability). Since there is only one response
  # variable, with only two elements, set intercept as 0.
  clin.effect[ncarr] ~ dnorm(0, 0.001)                # params = 1
  clin.effect[nclin] <- -clin.effect[ncarr]

  # equivalent to clin.effect + wt.effect:
  gp.effect ~ dnorm(clin.effect[gp_clinical], tau.wt) # params = 2
  v.effect ~ dnorm(clin.effect[v_clinical], tau.wt)   # params = 2
  o.effect ~ dnorm(clin.effect[o_clinical], tau.wt)   # params = 2

  # Prior distribution for wt.effect (log-odds for each ward type).
  # Sample different clin.effect from normal distribution for each clinical
  # class and convert to a probability).
  for (g in genders)
  {

    for (c in c(ncarr, nclin))
    {
      for (wt in hosp_wardtypes)
      {
        wtcg.effect[wt,c,g] <- clin.effect[c] + wt.effect[wt] + gender.effect[g]
      }

      for (w in hosp_wards)
      {
        cwtg_w[w,c,g] ~ dnorm(wtcg.effect[ward_type[w],c,g], tau.w)
        logit(cwtg_w.prob[w,c,g]) <- cwtg_w[w,c,g]      # params = 119*2*2
      }
    }

    # equivalent to cwtg_w.prob
    logit(gpg.prob[g]) <- gp.effect + gender.effect[g]
    logit(vg.prob[g]) <- v.effect + gender.effect[g]
    logit(og.prob[g]) <- o.effect + gender.effect[g]
  }

  # ------------------------

  # Prior value for intercept
  intercept ~ dnorm(0, 0.001)                         # params = 1

  # Prior values for precision
  tau.wt ~ dgamma(0.001, 0.001)                       # params = 1
  tau.w ~ dgamma(0.001, 0.001)                        # params = 1

  # Convert precisions to sd
  sd.wt <- 1/sqrt(tau.wt)
  sd.w <- 1/sqrt(tau.w)

  # Calculate odds
  c_diff <- clin.effect[ncarr] - clin.effect[nclin]
  odds_c <- exp(c_diff)
  g_diff <- gender.effect[male] - gender.effect[female]
  odds_g <- exp(g_diff)

  #monitor# full.pd, dic, deviance, gpg.prob, vg.prob, og.prob, odds_c, c_diff, g_diff, odds_g, sd.wt, sd.w
}

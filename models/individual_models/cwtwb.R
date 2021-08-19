model {
  # Define likelihood model for data:
  # Carbapenem resistance in hospital (gp, volunteer, and outpatient) samples
  # is Bernoulli distributed with probability wcb.prob (gpb.prob, vb.prob,
  # and ob.prob)
  for (p in 1:N_patients)
  {
    h_resist[p] ~ dbern(wcb_prob[ward[h_sample_GUID[p]],
                                 clinical[sample_type[h_sample_GUID[p]]],
				 h_bacteria[p]])
  }

  for (gp in 1:N_gp)
  {
    gp_resist[gp] ~ dbern(gpb_prob[h_bacteria[gp]])
  }

  for (v in 1:N_volunteers)
  {
    v_resist[v] ~ dbern(vb_prob[h_bacteria[v]])
  }

  for (o in 1:N_outpatients)
  {
    o_resist[o] ~ dbern(ob_prob[h_bacteria[o]])
  }

  # ------------------------

  # Carriage or Clinical
  mu_c[ncarr] ~ dnorm(0, tau)
  mu_c[nclin] <- -mu_c[ncarr]
  for (c in c(ncarr, nclin))
  {
    # A value for each hospital
    for (wt in hosp_wardtypes)
    {
      wtc_value[wt,c] ~ dnorm(mu_c[c], tau2)
      logit(wtc_prob[wt,c]) <- wtc_value[wt,c]
    }

    # A value for each ward
    for (w in hosp_wards)
    {
      wc_value[w,c] ~ dnorm(wtc_value[ward_type[w],c], tau3)
      for (b in bact_species)
      {
        logit(wcb_prob[w,c,b]) <- wc_value[w,c] + mu_b[b]
      }
    }
  }

  for (b in bact_species)
  {
    mu_b[b] ~ dnorm(mu, tau.b)
    logit(gpb_prob[b]) <- gp_value + mu_b[b]
    logit(vb_prob[b]) <- v_value + mu_b[b]
    logit(ob_prob[b]) <- o_value + mu_b[b]
  }

  # One for GP samples
  gp_value ~ dnorm(mu_c[gp_clinical], tau2)
  # One for volunteers
  v_value ~ dnorm(mu_c[v_clinical], tau2)
  # And one for outpatients
  o_value ~ dnorm(mu_c[o_clinical], tau2)

  # ------------------------

  tau ~ dgamma(0.001, 0.001)
  sd <- 1/sqrt(tau)
  tau2 ~ dgamma(0.001, 0.001)
  sd2 <- 1/sqrt(tau2)
  tau3 ~ dgamma(0.001, 0.001)
  sd3 <- 1/sqrt(tau3)
  tau.b ~ dgamma(0.001, 0.001)
  sd.b <- 1/sqrt(tau.b)

  mu ~ dnorm(0, 0.001)
  c_diff <- mu_c[ncarr] - mu_c[nclin]
  odds_c <- exp(c_diff)

  #monitor# full.pd, dic, deviance, sd, sd2, sd3, mu_c, mu_b, mu, odds_c, c_diff
}

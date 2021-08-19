model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~ dbern(asmagghwtw.prob[
        a,
        sample_month[h_sample_GUID[p]],
        age_group[h_sample_GUID[p]],
        gender[h_sample_GUID[p]],
        hospital[ward[h_sample_GUID[p]]],
        ward_type[ward[h_sample_GUID[p]]],
        ward[h_sample_GUID[p]]])
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~ dbern(asmagg.gpw.prob[
        a,
        sample_month[gp_sample_GUID[gp]],
        age_group[gp_sample_GUID[gp]],
        gender[gp_sample_GUID[gp]]])
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~ dbern(asmagg.vw.prob[
        a,
        sample_month[v_sample_GUID[v]],
        age_group[v_sample_GUID[v]],
        gender[v_sample_GUID[v]]])
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(asmagg.ow.prob[
        a,
        sample_month[o_sample_GUID[o]],
        age_group[o_sample_GUID[o]],
        gender[o_sample_GUID[o]]])
    }
  }

  # ------------------------

  # Define the priors:
  for(a in 1:antibiotic_classes)
  {
    antibiotic.class.effect[a] ~ dnorm(intercept, tau.class)
    logit(a.prob[a]) <- antibiotic.class.effect[a]
  }

  for(m in 1:N_sample_month)
  {
    samplemonth.effect[m] ~ dnorm(0, tau.samplemonth)
  }

  for (g in 1:N_age_group)
  {
    agegroup.effect[g] ~ dnorm(0, tau.agegroup)
    logit(agegroup.prob[g]) <- agegroup.effect[g]
  }

  gender.effect[female] ~ dnorm(0, 0.001)
  gender.effect[male] <- -gender.effect[female]

  for (h in 1:N_hosp)
  {
    h.effect[h] ~ dnorm(0, tau.hosp)

    for (wt in hosp_wardtypes)
    {
      hwt.effect[h,wt] ~ dnorm(h.effect[h], tau.wt)

      for (w in hosp_wards)
      {
        hwtw.effect[h,wt,w] ~ dnorm(hwt.effect[h,wt], tau.ward)
      }
    }
  }

  # equivalent to h.effect
  nh.effect ~ dnorm(0, tau.hosp)

  # equivalent to hwt.effect
  gp.effect ~ dnorm(nh.effect, tau.wt)
  v.effect ~ dnorm(nh.effect, tau.wt)
  o.effect ~ dnorm(nh.effect, tau.wt)

  # equivalent to hwtw.effect
  gpw.effect ~ dnorm(gp.effect, tau.ward)
  vw.effect ~ dnorm(v.effect, tau.ward)
  ow.effect ~ dnorm(o.effect, tau.ward)

  # ------------------------

  for(a in 1:antibiotic_classes)
  {
    for(m in 1:N_sample_month)
    {
      for (r in 1:N_age_group)
      {
        for (g in genders)
        {
          asmagg.effect[a,m,r,g] <-  antibiotic.class.effect[a] +
            samplemonth.effect[m] + agegroup.effect[r] + gender.effect[g]
          for (h in 1:N_hosp)
          {
            for (wt in hosp_wardtypes)
            {
              for (w in hosp_wards)
              {
                logit(asmagghwtw.prob[a,m,r,g,h,wt,w]) <-
                  asmagg.effect[a,m,r,g] +
                  hwtw.effect[h,wt,w]
              }
            }
          }

          logit(asmagg.gpw.prob[a,m,r,g]) <- asmagg.effect[a,m,r,g] + gpw.effect
          logit(asmagg.vw.prob[a,m,r,g]) <- asmagg.effect[a,m,r,g] + vw.effect
          logit(asmagg.ow.prob[a,m,r,g]) <- asmagg.effect[a,m,r,g] + ow.effect
        }
      }
    }
  }

  # ------------------------

  # Prior value for intercept (log-odds of the average resistance in all samples)
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)
  tau.samplemonth ~ dgamma(0.001, 0.001)
  tau.agegroup ~ dgamma(0.001, 0.001)
  tau.hosp ~ dgamma(0.001, 0.001)
  tau.wt ~ dgamma(0.001, 0.001)
  tau.ward ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.samplemonth <- sqrt(1/tau.samplemonth)
  sd.agegroup <- sqrt(1/tau.agegroup)
  sd.hosp <- sqrt(1/tau.hosp)
  sd.wt <- sqrt(1/tau.wt)
  sd.ward <- sqrt(1/tau.ward)

  #monitor# full.pd, dic, deviance, a.prob, intercept, sd.class, sd.samplemonth, sd.agegroup, sd.hosp, sd.wt, sd.ward
}

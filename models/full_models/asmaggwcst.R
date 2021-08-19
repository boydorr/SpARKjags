model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~
        dbern(asmaggwst.prob[a,
                             sample_month[h_sample_GUID[p]],
                             age_group[h_sample_GUID[p]],
                             gender[h_sample_GUID[p]],
                             ward[h_sample_GUID[p]],
                             clinical[sample_type[h_sample_GUID[p]]],
                             sample_type[h_sample_GUID[p]]])
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~ dbern(gp.prob[a,
                                              sample_month[gp_sample_GUID[gp]],
                                              age_group[gp_sample_GUID[gp]],
                                              gender[gp_sample_GUID[gp]],
                                              sample_type[gp_sample_GUID[gp]]])
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~ dbern(v.prob[a,
                                           sample_month[v_sample_GUID[v]],
                                           age_group[v_sample_GUID[v]],
                                           gender[v_sample_GUID[v]],
                                           sample_type[v_sample_GUID[v]]])
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(o.prob[a,
                                           sample_month[o_sample_GUID[o]],
                                           age_group[o_sample_GUID[o]],
                                           gender[o_sample_GUID[o]],
                                           sample_type[o_sample_GUID[o]]])
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

  for (w in hosp_wards)
  {
    ward.effect[w] ~ dnorm(0, tau.ward)
  }

  # equivalent to clin.effect + ward.effect
  gpcw.effect ~ dnorm(gp.effect, tau.ward)
  vcw.effect ~ dnorm(gp.effect, tau.ward)
  ocw.effect ~ dnorm(gp.effect, tau.ward)

  clin.effect[ncarr] ~ dnorm(0, 0.001)
  clin.effect[nclin] <- -clin.effect[ncarr]
  gp.effect <- clin.effect[gp_clinical]
  o.effect  <- clin.effect[o_clinical]
  v.effect  <- clin.effect[v_clinical]

  for (st in sampletypes)
  {
    st.effect[st] ~ dnorm(0, tau.st)
  }

  # ------------------------

  for(a in 1:antibiotic_classes)
  {
    for(m in 1:N_sample_month)
    {
      for (r in 1:N_age_group)
      {
        for (g in genders)
        {
          for (st in sampletypes)
          {
            asmaggst.effect[a,m,r,g,st] <- antibiotic.class.effect[a] +
              samplemonth.effect[m] + agegroup.effect[r] +
              gender.effect[g] + st.effect[st]

            for(w in hosp_wards)
            {
              for (c in c(ncarr, nclin))
              {
                logit(asmaggwst.prob[a,m,r,g,w,c,st]) <-
                  asmaggst.effect[a,m,r,g,st] + ward.effect[w] + clin.effect[c]
              }
            }

            logit(gp.prob[a,m,r,g,st]) <- asmaggst.effect[a,m,r,g,st] +
              gpcw.effect
            logit(v.prob[a,m,r,g,st]) <- asmaggst.effect[a,m,r,g,st] +
              vcw.effect
            logit(o.prob[a,m,r,g,st]) <- asmaggst.effect[a,m,r,g,st] +
              ocw.effect
          }
        }
      }
    }
  }

  # Prior value for intercept (log-odds of the average resistance in all samples)
  intercept ~ dnorm(0, 0.001)

  # Prior values for precision
  tau.class ~ dgamma(0.001, 0.001)
  tau.samplemonth ~ dgamma(0.001, 0.001)
  tau.agegroup ~ dgamma(0.001, 0.001)
  tau.ward ~ dgamma(0.001, 0.001)
  tau.st ~ dgamma(0.001, 0.001)

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.samplemonth <- sqrt(1/tau.samplemonth)
  sd.agegroup <- sqrt(1/tau.agegroup)
  sd.ward <- sqrt(1/tau.ward)
  sd.st <- sqrt(1/tau.st)

  #monitor# full.pd, dic, deviance, a.prob, intercept, sd.class, sd.samplemonth, sd.agegroup, sd.ward, sd.st
}

model {
  # Define likelihood model for data:
  for (p in 1:N_patients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[h_GUID[p],a] ~ dbern(aagg.prob[a,
                                              sample_month[h_sample_GUID[p]],
                                              gender[h_sample_GUID[p]],
                                              age_group[h_sample_GUID[p]]])
    }
  }

  for (gp in 1:N_gp)
  {
    for(a in 1:antibiotic_classes)
    {
      response[gp_GUID[gp],a] ~ dbern(aagg.prob[a,
                                                sample_month[gp_sample_GUID[gp]],
                                                gender[gp_sample_GUID[gp]],
                                                age_group[gp_sample_GUID[gp]]])
    }
  }

  for (v in 1:N_volunteers)
  {
    for(a in 1:antibiotic_classes)
    {
      response[v_GUID[v],a] ~ dbern(aagg.prob[a,
                                              sample_month[v_sample_GUID[v]],
                                              gender[v_sample_GUID[v]],
                                              age_group[v_sample_GUID[v]]])
    }
  }

  for (o in 1:N_outpatients)
  {
    for(a in 1:antibiotic_classes)
    {
      response[o_GUID[o],a] ~ dbern(aagg.prob[a,
                                              sample_month[o_sample_GUID[o]],
                                              gender[o_sample_GUID[o]],
                                              age_group[o_sample_GUID[o]]])
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

  for(a in 1:antibiotic_classes)
  {
    for(m in 1:N_sample_month)
    {
    for (g in genders)
    {
      for (r in 1:N_age_group)
      {
        logit(aagg.prob[a,m,g,r]) <- antibiotic.class.effect[a] +
          samplemonth.effect[m] + agegroup.effect[r] + gender.effect[g]
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

  # Convert precisions to sd
  sd.class <- sqrt(1/tau.class)
  sd.samplemonth <- sqrt(1/tau.samplemonth)
  sd.agegroup <- sqrt(1/tau.agegroup)

  #monitor# full.pd, dic, deviance, a.prob, intercept, sd.agegroup, sd.samplemonth, sd.class
}

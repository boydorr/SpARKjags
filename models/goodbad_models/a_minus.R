#' $$y_{ij} ~ Bern(a_{g})$$
#' where $y_{ij}$ is a binary variable denoting resistance of the $i^{th}$
#' isolate to the $a^{th}$ antibiotic class and $a_g$ is

model {
  # Define likelihood model for data, $P(\theta|data)$ ----------------------

  # Stochastic component of likelihood function, linking the response
  # variable to `a.prob`, given Bernoulli distributed sampling error
  for (p in 1:N_patients) {
    # Probability of belonging to the bad group
    bad.p[p] ~ dbern(prob.of.bad.hosp)
    # Index the bad group in a.prob
    index.bad.p[p] <- bad.p[p] + 1

    for(a in 1:antibiotic_classes) {
      # Response is different for each antibiotic and depending on
      # which pop it's from
      response[h_GUID[p],a] ~ dbern(a.prob[a, index.bad.p[p]])
    }
  }

  for (gp in 1:N_gp) {
    bad.gp[gp] ~ dbern(prob.of.bad.gp)
    index.bad.gp[gp] <- bad.gp[gp] + 1

    for(a in 1:antibiotic_classes) {
      response[gp_GUID[gp],a] ~ dbern(a.prob[a, index.bad.gp[gp]])
    }
  }

  for (v in 1:N_volunteers) {
    bad.v[v] ~ dbern(prob.of.bad.vol)
    index.bad.v[v] <- bad.v[v] + 1

    for(a in 1:antibiotic_classes) {
      response[v_GUID[v],a] ~ dbern(a.prob[a, index.bad.v[v]])
    }
  }

  for (o in 1:N_outpatients) {
    bad.o[o] ~ dbern(prob.of.bad.out)
    index.bad.o[o] <- bad.o[o] + 1

    for(a in 1:antibiotic_classes) {
      response[o_GUID[o],a] ~ dbern(a.prob[a, index.bad.o[o]])
    }
  }

  # A component to track `a.prob` predicted by the model
  for(a in 1:antibiotic_classes) {
    # probability of being resistant in the good group (less resistances)
    antibiotic.class.effect[a] ~ dnorm(intercept, tau.class)
    logit(a.prob[a,1]) <- antibiotic.class.effect[a]

    # probability of being resistant in the bad group (many resistances) will
    # always be higher than the good group
    logit(a.prob[a,2]) <- antibiotic.class.effect[a] + diff
  }

  # Define the priors, $P(\theta)$ ------------------------------------------

  # Prior value for intercept (make the precision small to emphasize the lack
  # of prior information)
  intercept ~ dnorm(0, 0.001)

  # Difference between distribution means of the "good" and "bad" effect
  diff ~ dgamma(0.001, 0.001)
  intercept.plus <- intercept + diff

  # Probability of being in the bad group
  prob.of.bad.hosp ~ dbeta(1, 1)
  prob.of.bad.gp ~ dbeta(1, 1)
  prob.of.bad.vol ~ dbeta(1, 1)
  prob.of.bad.out ~ dbeta(1, 1)

  # Estimate a measure of variation (precision) for the sampling error
  # distribution
  tau.class ~ dgamma(0.001, 0.001)

  # Convert precisions to SD ------------------------------------------------

  sd.class <- sqrt(1/tau.class)

  #monitor# full.pd, dic, deviance, a.prob, prob.of.bad.hosp, prob.of.bad.gp, prob.of.bad.vol, prob.of.bad.out, bad.p, bad.gp, bad.v, bad.o, intercept, sd.class
}

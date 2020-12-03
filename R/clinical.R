clinical <- function(data, what) {

  clinical <- data$lookup$clinical %>%
    dplyr::filter(clinical == "yes")
  index <- clinical$index

  # Return the index associated with clinical
  if(what == "index") {
    return(index)

    # Return whether gp is clinical
  }else if(what == "gp") {
    gp <- data$jags$gp_clinical
    return(gp == index)

    # Return whether outpatient is clinical
  }else if(what == "out") {
    out <- data$jags$o_clinical
    return(out == index)

    # Return whether volunteer is clinical
  }else if(what == "vol") {
    vol <- data$jags$v_clinical
    return(vol == index)
  }else {
    stop("unknown")
  }
}

carriage <- function(data, what) {

    carriage <- data$lookup$clinical %>%
      dplyr::filter(clinical == "no")
    index <- carriage$index

    # Return the index associated with clinical
    if(what == "index") {
      return(index)

      # Return whether gp is clinical
    }else if(what == "gp") {
      gp <- data$jags$gp_clinical
      return(gp == index)

      # Return whether outpatient is clinical
    }else if(what == "out") {
      out <- data$jags$o_clinical
      return(out == index)

      # Return whether volunteer is clinical
    }else if(what == "vol") {
      vol <- data$jags$v_clinical
      return(vol == index)
    }else {
      stop("unknown")
    }
}

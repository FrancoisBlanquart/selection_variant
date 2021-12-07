Rtor <- function(R, mu, SD){
  V <- SD^2
  return(
    (mu/V) * (R^(V/mu^2) - 1) 
  )
}
rtoR <- function(r, mu, SD){
  V <- SD^2
  return(
    (1 + (V/mu) * r)^(mu^2/V)
  )
}
get_gamma_function <- function(mu, sd){ # get parameters of gamma distribution with mean mu and sd sd
  return(list(shape = mu^2/sd^2, scale = sd^2/mu))
}
# get means of the var-cov matrix describing the distribution (rH, rV) (the means are the true rH and rV for each timepoint)
getmeans <- function(pars, Rpars, weekly = F){
  
  if(weekly) weekly_factor <- 7 else weekly_factor <- 1
  # stopifnot(length(Rpars)==ndays-1) # R of wild-type each day
  s1 <- pars[1]
  s2 <- pars[2]
  s3 <- pars[3]
  rH <- Rtor(Rpars, MUW/weekly_factor, SDW/weekly_factor)
  rV <- Rtor(Rpars * (1+s1), MUW*(1+s2)/weekly_factor, SDW * (1+s3)/weekly_factor)
  return(c(rbind(rH, rV)))
  
}
getvar <- function(i, tab){ # variance of rH and rV
  nn <- nrow(tab)
  stopifnot(i>=2 & i<=nn)
  # i: time index
  I0 <- tab$Psmoothed[i-1] # number of cases
  I1 <- tab$Psmoothed[i]
  N0 <- tab$Ntrue[i-1] # number of cases assayed for VOC
  N1 <- tab$Ntrue[i]
  p0 <- tab$fVOC[i-1] # frequency of variant
  p1 <- tab$fVOC[i]
  return(
    list(
      vH = 1/I1 + 1/I0 + p0/(N0*(1-p0)) + p1/(N1*(1-p1)),
      vM = 1/I1 + 1/I0 + (1-p0)/(p0*N0) + (1-p1)/(p1*N1)
    )
  )
}
getcov <- function(i, tab){ # covariance of RH and rM at same time-point
  nn <- nrow(tab)
  stopifnot(i>=2 & i<=nn)
  # i: time index
  I0 <- tab$Psmoothed[i-1] # number of cases
  I1 <- tab$Psmoothed[i]
  N0 <- tab$Ntrue[i-1] # number of cases assayed for VOC
  N1 <- tab$Ntrue[i]
  return(
    1/I1 + 1/I0 - 1/N0 - 1/N1
  )
}
getspatialvar <- function(i, tab){ # variance of rH and rV
  nn <- nrow(tab)
  stopifnot(i>=2 & i<=nn)
  # i: country index
  I0 <- tab$Psmoothed0[i] # number of cases
  I1 <- tab$Psmoothed1[i]
  N0 <- tab$Ntrue0[i] # number of cases assayed for VOC
  N1 <- tab$Ntrue1[i]
  p0 <- tab$fVOC0[i] # frequency of variant
  p1 <- tab$fVOC1[i]
  return(
    list(
      vH = 1/I1 + 1/I0 + p0/(N0*(1-p0)) + p1/(N1*(1-p1)),
      vM = 1/I1 + 1/I0 + (1-p0)/(p0*N0) + (1-p1)/(p1*N1)
    )
  )
}
getspatialcov <- function(i, tab){ # covariance of RH and rM at same time-point
  nn <- nrow(tab)
  stopifnot(i>=2 & i<=nn)
  # i: country index
  I0 <- tab$Psmoothed0[i] # number of cases
  I1 <- tab$Psmoothed1[i]
  N0 <- tab$Ntrue0[i] # number of cases assayed for VOC
  N1 <- tab$Ntrue1[i]
  return(
    1/I1 + 1/I0 - 1/N0 - 1/N1
  )
}
getcrosscovH <- function(i, tab){ # covariance between rH and rH of successive timepoints
  nn <- nrow(tab)
  # covariance between r of times i andgm i+1
  stopifnot(i>=2 & i<=nn-1)
  # i: time index
  I1 <- tab$Psmoothed[i]
  N1 <- tab$Ntrue[i]
  p1 <- tab$fVOC[i]
  return(
    -1/I1 - p1/(N1 * (1-p1))
  )
}
getcrosscovV <- function(i, tab){ # covariance between rV and rV of successive timepoints
  nn <- nrow(tab)
  # covariance between r of times i and i+1
  stopifnot(i>=2 & i<=nn-1)
  # i: time index
  I1 <- tab$Psmoothed[i]
  N1 <- tab$Ntrue[i]
  p1 <- tab$fVOC[i]
  return(
    -1/I1 - (1-p1)/(N1 * p1)
  )
}
getcrosscovHV <- function(i, tab){ # covariance between rH and rV of successive timepoints
  nn <- nrow(tab)
  # covariance between r of times i and i+1
  stopifnot(i>=2 & i<=nn-1)
  # i: time index
  I1 <- tab$Psmoothed[i]
  N1 <- tab$Ntrue[i]
  return(
    -1/I1 + 1/N1
  )
}
create_cov_mat <- function(tab){
  nn <- nrow(tab)
  stopifnot(all(c("Psmoothed", "Ntrue", "fVOC") %in% names(tab)))
  
  varcov <- matrix(0, nrow = 2*(nn-1), ncol = 2*(nn-1))
  diag(varcov) <- 1
  
  times <- c(sapply(2:nn, rep, 2))
  vars <- rep(c("H", "V"), nn-1)
  
  # NOW FILL IN THE VARCOV MATRIX:
  # 1) diagonal
  diag(varcov) <- unlist(c(sapply(2:nn, getvar, tab = tab)))
  
  # 2) rHt and rVt covariances
  for(i in 2:nn){
    idx <- which(times==i)
    stopifnot(length(idx)==2)
    stopifnot(vars[idx][1] == "H" & vars[idx][2] == "V")
    varcov[idx[1], idx[2]] <- getcov(i, tab = tab)
    varcov[idx[2], idx[1]] <- getcov(i, tab = tab)
  }
  # 2) rHt and rHt+1 covariances
  for(i in 2:(nn-1)){
    idx0 <- which(times==i & vars == "H")
    idx1 <- which(times==i+1 & vars == "H")
    stopifnot(length(idx0)==1); stopifnot(length(idx1) == 1)
    varcov[idx0, idx1] <- getcrosscovH(i, tab = tab)
    varcov[idx1, idx0] <- getcrosscovH(i, tab = tab)
  }
  # 3) rMt and rMt+1 covariances
  for(i in 2:(nn-1)){
    idx0 <- which(times==i & vars == "V")
    idx1 <- which(times==i+1 & vars == "V")
    stopifnot(length(idx0)==1); stopifnot(length(idx1) == 1)
    varcov[idx0, idx1] <- getcrosscovV(i, tab = tab)
    varcov[idx1, idx0] <- getcrosscovV(i, tab = tab)
  }
  # 4) fill in the rMt and rHt+1 covariances
  for(i in 2:(nn-1)){
    # rMt and rHt+1
    idx0 <- which(times==i & vars == "V")
    idx1 <- which(times==i+1 & vars == "H")
    stopifnot(length(idx0)==1); stopifnot(length(idx1) == 1)
    varcov[idx0, idx1] <- getcrosscovHV(i, tab = tab)
    varcov[idx1, idx0] <- getcrosscovHV(i, tab = tab)
    # rHt and rMt+1
    idx0 <- which(times==i+1 & vars == "V")
    idx1 <- which(times==i & vars == "H")
    stopifnot(length(idx0)==1); stopifnot(length(idx1) == 1)
    varcov[idx0, idx1] <- getcrosscovHV(i, tab = tab)
    varcov[idx1, idx0] <- getcrosscovHV(i, tab = tab)
    
  }
  return(varcov)
}
create_spatial_cov_mat <- function(tab){
  
  nn <- nrow(tab)
  stopifnot(all(c("Psmoothed0", "Psmoothed1", "Ntrue0", "Ntrue1", "fVOC0", "fVOC1") %in% names(tab)))
  
  varcov <- matrix(0, nrow = 2*(nn-1), ncol = 2*(nn-1))
  diag(varcov) <- 1
  
  countries <- c(sapply(2:nn, rep, 2))
  vars <- rep(c("H", "V"), nn-1)
  
  # NOW FILL IN THE VARCOV MATRIX:
  # 1) diagonal
  diag(varcov) <- unlist(c(sapply(2:nn, getspatialvar, tab = tab)))
  
  # 2) rHt and rVt covariances
  for(i in 2:nn){
    idx <- which(countries==i)
    stopifnot(length(idx)==2)
    stopifnot(vars[idx][1] == "H" & vars[idx][2] == "V")
    varcov[idx[1], idx[2]] <- getspatialcov(i, tab = tab)
    varcov[idx[2], idx[1]] <- getspatialcov(i, tab = tab)
  }
  return(varcov)
}
mylik <- function(allpars, mydata, vars, weekly, return_full = F){
  
  # parameters are s1, s2, s3 fixed to 0
  # then all rW (length ndays-1)
  
  #stopifnot(length(allpars) == npar1)
  means <- getmeans(c(allpars[1], allpars[2], allpars[3]), Rpars = allpars[4:length(allpars)], weekly = weekly)
  
  #print(means)
  # if(weekly){
  #   vars <- varcov_weekly
  # } else {
  #   vars <- varcov
  # }
  #print(means)
  
  stopifnot(length(means) == nrow(vars) & length(means) == ncol(vars)) # check indeed the mean is same length as the dimension of variance - covariance matrix
  
  if(return_full){
    return(list(
      lik = mnormt::dmnorm(x = mydata, mean = means, varcov = vars, log = T),
      means = means
    ))
  } else {
    mnormt::dmnorm(x = mydata, mean = means, varcov = vars, log = T) 
  }
}

fastneglik <- function(allpars, mydata, inverse_vars, logdetvars, fixed_s3, weekly){
  n <- (length(allpars) - 2) * 2 # dimension of the matrix = number of R times 2
  means <- getmeans(c(allpars[1], allpars[2], fixed_s3), Rpars = allpars[3:length(allpars)], weekly = weekly)
  stopifnot(n==nrow(inverse_vars))
  return(
    -as.numeric(unname(
      unlist(
      -n/2 * log(2*pi) - (1/2) * logdetvars - (1/2) * matrix(mydata - means, nrow = 1) %*% inverse_vars %*% matrix(mydata - means, nrow = n)
    )))
  )
}
gradient <- function(allpars, mydata, inverse_vars, logdetvars, weekly){
  
}

myneglik <- function(allpars, mydata, vars, fixed_s3, weekly, return_full) - mylik(c(allpars[1:2], fixed_s3, allpars[3:length(allpars)]), mydata, vars, weekly, return_full)

mylik2 <- function(mypars, mydata, myinitRW, vars, fixed_s3, weekly) mylik(allpars = c(mypars, fixed_s3, myinitRW), mydata, vars, weekly) # effect on variance set to 0, RW set to as estimated from data
myneglik2 <- function(mypars, mydata, myinitRW, vars, fixed_s3, weekly) - mylik2(mypars, mydata, myinitRW, vars, fixed_s3, weekly)

fastneglik2 <- function(mypars, mydata, myinitRW, inverse_vars, logdetvars, fixed_s3, weekly) fastneglik(allpars = c(mypars, myinitRW), mydata, inverse_vars, logdetvars, fixed_s3, weekly)

deriv_matrix <- function(mypars, myinitRW){
  # matrix of partial derivatives of means wrt parameters, used to compute gradients
  npar <- length(allpars)
  n <- (length(allpars) - 2) * 2 # number of rows in matrix
  m <- matrix(NA, nrow = n, ncol = npar)
  print("this function is unfinished!")
  stop()
}

create_BMautocov_mat <- function(sigma, max_time){
  BMautocov_mat <- matrix(NA, nrow = max_time, ncol = max_time)
  for(i in 1:max_time){
    for(j in 1:max_time){
      if(i == j) BMautocov_mat[i, j] <- sigma^2 * i
      if(i != j) BMautocov_mat[i, j] <- sigma^2 * min(i, j)
    }
  }
  return(BMautocov_mat)
}
prior_RW <- function(allpars){
  stop("must adapt depending on daily or weekly data")
  # Attempt at a prior (what variance, what temporal autocorrelation to choose?)
  # Brownian motion: covariance between two values B_s and B_t is the variance of the process x min(s,t)
  Rpars <- allpars[4:(ndays+2)] # R parameters
  mnormt::dmnorm(x = Rpars, mean = rep(mean_RW, ndays-1), varcov = BMautocov_mat, log = T) # imposed temporal autocovariance structure
}

# generate simulated data in same format as our data
generate_sim_data <- function(mysim, freq_sample_size, weekly){
  
  mysim_sub_b <- data.frame(
    day = 0:(max_time-1),
    data = "SIM",
    Psmoothed = mysim$all_detected_incidenceW+mysim$all_detected_incidenceM,
    Psmoothed_H = mysim$all_detected_incidenceW,
    Psmoothed_V = mysim$all_detected_incidenceM,
    N = freq_sample_size # sample size for frequency estimation
  )
  mysim_sub_b$week <- floor(mysim_sub_b$day/7) + 1
  if(weekly){ # aggregate by week if necessary
    mysim_sub_b <- ddply(mysim_sub_b, .(week), summarise, week = unique(week), data =  unique(data), Psmoothed = sum(Psmoothed), Psmoothed_H = sum(Psmoothed_H), Psmoothed_V = sum(Psmoothed_V), N = sum(N))
  }
  mysim_sub_b$fVOC <- mysim_sub_b$Psmoothed_V / (mysim_sub_b$Psmoothed_H + mysim_sub_b$Psmoothed_V)
  
  nn <- nrow(mysim_sub_b) # number of rows now depends on weekly or not
  if(weekly) n_to_keep <- nn-burnin else n_to_keep <- nn-burnin # how many points to keep
    
  # generate number of cases and frequencies with error:
  mysim_sub_b$Psmoothed_error <- rpois(n = nn, lambda = mysim_sub_b$Psmoothed) # Poisson cases
  mysim_sub_b$fVOC_error <- rbinom(n = rep(1, nn), size = mysim_sub_b$N, prob = mysim_sub_b$fVOC) / mysim_sub_b$N

  # generate exact growth rates
  #mysim_sub_b$Psmoothed_V <- mysim_sub_b$Psmoothed * mysim_sub_b$fVOC
  #mysim_sub_b$Psmoothed_H <- mysim_sub_b$Psmoothed * (1-mysim_sub_b$fVOC)
  mysim_sub_b$rH <- NA
  mysim_sub_b$rV <- NA
  mysim_sub_b[2:nn, "rH"] <- log(mysim_sub_b[2:nn, "Psmoothed_H"] / mysim_sub_b[1:(nn-1), "Psmoothed_H"])
  mysim_sub_b[2:nn, "rV"] <- log(mysim_sub_b[2:nn, "Psmoothed_V"] / mysim_sub_b[1:(nn-1), "Psmoothed_V"])

  # generate growth rates with error
  mysim_sub_b$Psmoothed_V_error <- mysim_sub_b$Psmoothed_error * mysim_sub_b$fVOC_error
  mysim_sub_b$Psmoothed_H_error <- mysim_sub_b$Psmoothed_error * (1-mysim_sub_b$fVOC_error)
  mysim_sub_b$rH_error <- NA
  mysim_sub_b$rV_error <- NA
  mysim_sub_b[2:nn, "rH_error"] <- log(mysim_sub_b[2:nn, "Psmoothed_H_error"] / mysim_sub_b[1:(nn-1), "Psmoothed_H_error"])
  mysim_sub_b[2:nn, "rV_error"] <- log(mysim_sub_b[2:nn, "Psmoothed_V_error"] / mysim_sub_b[1:(nn-1), "Psmoothed_V_error"])
  
  # get rid of initial (transient) effects:
  mysim_sub_b <- mysim_sub_b[(nn-n_to_keep+1):nn, ]
  
  return(mysim_sub_b)
}

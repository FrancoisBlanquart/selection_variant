mylogit <- function(p) log(p/(1-p))

myinvlogit <- function(x) exp(x) / (1 + exp(x))

myintegrate <- function(f, tmin, tmax, deltat, param, ...){ # my function to simply integrate a function with parameters "param"
  if(tmin == tmax) return(0)
  tmp <- sapply(seq(tmin, tmax, deltat), f, param = param, ...)
  if(any(is.infinite(tmp))) stop("function takes infinite values")
  n <- length(tmp)
  weights <- rep(2, n); weights[1] <- 1; weights[n] <- 1
  return(sum(0.5 * (tmax-tmin)/(n-1) * weights * tmp))
}

solve_rtoR <- function(param, r){
  if(param[1] < 1) stop("shape parameter must be greater than 1")
  f <- function(t, param, r) dgamma(t, shape = param[1], scale = param[2]) * exp(-r * t)
  return(
    1 /  myintegrate(f, tmin = 0, tmax = 40, delta = 0.0001, param = param, r = r)
  )
}

solve_Rtor <- function(param, R) uniroot(f = function(r) {solve_rtoR(param, r) - R}, interval = c(-1, 1))$root

# new version with explicit formula valid for gamma distribution
new_rtoR <- function(param, r) 1 / (1 + param[2] * r)^-param[1]
new_Rtor <- function(param, R) -1/param[2] + (1/param[2]) * (1/R)^(-1/param[1])

# test these functions
if(test_functions <- F){
  tmp <- c()
  for(i in 1:10){
    gtpar <- c(runif(1,1,10), runif(1, 0, 1))
    rrand <- runif(1, 0, 0.1)
    
    tmp <- rbind(tmp, c(solve_rtoR(param = gtpar, r = rrand),
                        new_rtoR(param = gtpar, r = rrand)))
  }
  print("range of error on the two functions from r to R:")
  print(range(abs((tmp[,2]-tmp[,1]))))
}

# weekly_to_r <- function(mf){
#   log(mf)/7
# }

getbetaparams <- function(mu, var) {
  # gets parameter for a beta distrib with mean mu and variance var
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
getlik_frequency <- function(NVOCdata, Ndata, daymin, daymax, simtable, myrho){ # get likelihood component TODO see if can make this faster
  
  # get simulated frequency of VOC in sim:
  idx <- which(simtable$dateday >= daymin  & simtable$dateday <= daymax)
  if(length(idx)<1) stop("data not present")
  
  fsim <- sum(simtable$PVOCsim[idx]) / sum(simtable$PTOTsim[idx]) 
  
  # probability to get NVOCdata VOC if true frequency is fsim
  return(
    - sum(dbetabinom(x = NVOCdata,  prob = fsim, rho = myrho, size = Ndata, log = T))
    # note to myself: the maximum likelihood estimator of beta binomial will not be fsim = NVOCdata/Ndata but higher
    #- sum(dbinom(x = NVOCdata,  prob = fsim, size = Ndata, log = T))
  )
}
regional_lik <- function(par, region, epitable, freqdatatable, return_full = FALSE){
  
  stopifnot(length(par)==8)
  # region specific parameters:
  p0 <- par[1]
  Rinit <- par[2]
  Rfinal <- par[3]
  k <- par[4] # > 0, not too large
  inflection_date <- par[5] # must be in minday:maxday
  
  # generic parameters
  transmission_advantage <- par[6]
  size_binomial <- par[7] # dispersion of negative binomial
  rho <- par[8] # dispersion of beta binomial
  
  RWT <- Rinit + (Rfinal - Rinit) / (1 + exp(- k * (minday:(maxday-1) - inflection_date)))
  RVOC <-  RWT * (1 + transmission_advantage)
  
  rWT <- new_Rtor(gtpar, RWT)
  rVOC <- new_Rtor(gtpar, RVOC)
  
  stopifnot(length(rWT) == n-1)
  stopifnot(length(RVOC) == n-1)
  
  epitable$PVOCsim <- NA
  epitable$PWTsim <- NA
  
  # initialise prevalence of VOC and WT using frequency parameter
  epitable[1, "PVOCsim"] <- p0 * epitable[1, "Psmoothed"] 
  epitable[1, "PWTsim"] <- (1-p0) * epitable[1, "Psmoothed"]
  
  # extrapolate with exp growth from minday+1 to maxday
  if(!return_full){
    epitable[2:n, "PWTsim"]  <- epitable[1 , "PWTsim"] * exp(cumsum(rWT))
    epitable[2:n, "PVOCsim"] <- epitable[1 , "PVOCsim"] * exp(cumsum(rVOC))
  } else {
    # add rs values equal to the last one:
    finalrWT <- rWT[n-1]
    finalrVOC <- rVOC[n-1]
    rWT <- c(rWT, rep(finalrWT, 90))
    rVOC <- c(rVOC, rep(finalrVOC, 90))
    epitable <- rbind.fill(epitable, expand.grid(reg2 = region, dateday = (maxday+1):(maxday+90)))
    epitable[2:nrow(epitable), "PWTsim"]  <- epitable[1 , "PWTsim"] * exp(cumsum(rWT))
    epitable[2:nrow(epitable), "PVOCsim"] <- epitable[1 , "PVOCsim"] * exp(cumsum(rVOC))
  }
  
  epitable$PTOTsim <- epitable$PVOCsim+epitable$PWTsim
  
  # component of the likelihood based on number of Psmoothed daily compared to inital prediction
  ll0 <- - sum(dnbinom(x = epitable$Psmoothed, size = size_binomial, mu = epitable$PTOTsim, log = T), na.rm = T)
  
  # LIKELIHOOD COMPONENT CORRESPONDING TO FREQUENCIES
  
  if(nrow(freqdatatable) > 0){ # if frequency data is available
    lik_freqs <- sum(
      apply(freqdatatable, 1, function(vec) getlik_frequency(NVOCdata = as.numeric(vec["NVOC"]), Ndata = as.numeric(vec["N"]), daymin = as.numeric(vec["daymin"]), daymax = as.numeric(vec["daymax"]), simtable = epitable, myrho = rho))
    )
  }
  
  if(return_full){
    return(list(ll = ll0+sum(lik_freqs), epitable = epitable))
  } else {
    return(ll0+sum(lik_freqs))
  }
}
optim.fun.repeated <- function(n.repeats, lik.fun, init.fun, optim.fun = "nmk", verbose = F, lower = NULL, upper = NULL, ...){
  
  # a function to MINIMISE properly the negative log-likelihood function
  # using repeated use of function optim
  
  generate.init.fun <- match.fun(init.fun)
  lik.fun <- match.fun(lik.fun)
  optim.fun.text <- optim.fun; optim.fun <- match.fun(optim.fun.text) # optimisation function, can be nmk/nmkb or optim
  
  all.opt <- list()
  for(ii in 1:n.repeats){
    init <- generate.init.fun()
    
    if(length(init)==1){ # one parameter function, use optimise
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use (0,1)
        interval <- c(0, 1)
      } else {
        interval <- c(lower[1], upper[1])
      }
      all.opt[[ii]] <- optimise(lik.fun, interval = interval, ...)
      names(all.opt[[ii]]) <- c("par", "value")
      all.opt[[ii]]$message <- "Successful convergence"
    } else { # otherwise use nmk
      
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use nmk
        if(optim.fun.text == "nmk") all.opt[[ii]] <- optim.fun(init, lik.fun, ...)
        if(optim.fun.text == "optim") all.opt[[ii]] <- optim.fun(init, lik.fun, method = "BFGS", ...)
        if(optim.fun.text != "nmk" & optim.fun.text != "optim") stop()
        flag <- F
      } else {
        all.opt[[ii]] <- NA
        while(is.na(all.opt[[ii]])){
          init <- generate.init.fun()
          if(optim.fun.text == "nmkb") all.opt[[ii]] <- optim.fun(init, lik.fun, lower = lower, upper = upper, ...)
          if(optim.fun.text == "optim") all.opt[[ii]] <- optim.fun(init, lik.fun, method = "BFGS", lower = lower, upper = upper, ...)
          if(optim.fun.text != "nmkb" & optim.fun.text != "optim") stop()
        }
      }
    }
    if(verbose)print(ii)
  }
  if(verbose) cat(sum(unlist(lapply(all.opt, function(x)x$message == "Successful convergence"))), "optimisations converge\n")
  all.ll <- unlist(lapply(all.opt, function(x)x$value))
  min(all.ll) -> minll.pow
  if(verbose) cat(length(which(all.ll < minll.pow + 1 & all.ll > minll.pow)), " optim within 1 log lik of minimum likelihood\n")
  output <- list(minll.pow, all.opt[[which.min(all.ll)]]$par)
  names(output) <- c("lik", "pars")
  return(output)
}
run_MCMC_Gibbs <- function(start.value, N.iter, proposal.sd, posterior.fun, npar, lower, upper, verbose = T, ...){
  # Gibbs sampling
  # log-normal proposal
  # code adapted from Henrik Salje, Cécile Tran Kiem  on https://zenodo.org/record/3813815#.YBx5MC1Q0Ut
  # explanations here for the log-normal proposal  https://umbertopicchini.wordpress.com/2017/12/18/tips-for-coding-a-metropolis-hastings-sampler/

  chain = array(dim = c(N.iter, npar))
  all.lik <- rep(NA, N.iter)
  acceptance <- matrix(0, nrow = N.iter, ncol = npar)
  chain[1,] = start.value
  all.lik[1] <- posterior.fun(start.value, ...)

  for (i in 2:N.iter){
    
    if(i / 1000 == floor(i / 1000) & verbose){
      print(i)
      acceptance_byparam <- colMeans(acceptance[1:i, ]) # assess acceptance rate
      cat("acceptance rate in last 1000 generations: ", acceptance_byparam, "\n")
      cat("likelihood: ", all.lik[i-1], "\n")
      #for(iparam in 1:npar) plot(chain[, iparam], type = "l", main = iparam)
    }
    
    # Getting parameters and posterior obtained at the previous iteration
    chain[i,]  <-  chain[i-1,]
    all.lik[i] <- all.lik[i-1]
    
    # Parameter updates
    for(iparam in 1:npar){
      
      old_param <- chain[i - 1, iparam]
      
      # sampling a new candidate in log-normal distribution with sd proposal.sd
      new_param <- old_param * exp(proposal.sd[iparam] * rnorm(1))
      
      if(new_param > upper[iparam] || new_param < lower[iparam]){ # if not in the limits, just keep old parameter
        chain[i, iparam] <- old_param
      } else { # if in the limits
        
        chain[i, iparam] <- as.numeric(new_param)
        
        # posterior function for this 
        newlik <- posterior.fun(chain[i,], ...)
        
        # acceptance of rejection
        log_acceptance_ratio <- newlik - all.lik[i] + log(new_param) - log(old_param) # the latter factor is a correction for the log-normal proposal
        
        if(log(runif(1)) < log_acceptance_ratio){
          all.lik[i] <- newlik
          acceptance[i, iparam] <- 1
        } else{
          chain[i, iparam] <- old_param
        }
      }
    }
  }
  return(list(chain = chain, all.lik = all.lik))
}
# testing the new MCMC function:
# fakedata <- rnorm(1000, mean = 1, sd = 3)
# mypost <- function(par) sum(dnorm(x = fakedata, mean = par[1], sd = par[2], log = T))
# mcmc <- run_MCMC_Gibbs(start.value = c(2, 1), N.iter = 10000, proposal.sd = c(0.1, 0.1), posterior.fun = mypost, npar = 2, lower=  c(0,0), upper = c(10, 10), verbose = T)
# plot(mcmc$all.lik)
# plot(mcmc$chain[,1])
# plot(mcmc$chain[,2])

run_MCMC_MH <- function(start.value, N.iter, proposal.cov, posterior.fun, npar, lower, upper, verbose = T, ...){ # same as run_MCMC but tuning the proposal distribution
  
  # start.value = starting values of parameters
  # N.iter = number of iterations
  # proposal.cov = covariance-matrix of the multivariate normal distribution used as proposal. Has dimension npar x npar. Usually starts with a diagonal matrix with chosen variances for each parameter
  # posterior.fun = posterior distribution (a function of the vector of parameters)
  # npar = number of parameters
  # lower = vector of lower bounds on parameters (length npar)
  # upper = vector of upper bounds on parameters (length npar)
  # ...  = extra arguments that need to be passed on to the posterior.fun (beyond the vector of parameters)
  
  # stop("verify algo; convergence problems?")
  # THIS IS THE LATEST AND MOST CORRECT VERSION OF THE MCMC ALGO
  # Metropolis-Hastings algorithm (requires symmetric proposal distribution -> multivariate Gaussian OK)
  init.proposal.cov <- proposal.cov
  factor.cov <- 10
  
  posterior.fun <- match.fun(posterior.fun)
  stopifnot(dim(proposal.cov) == c(npar, npar))
  
  chain = array(dim = c(N.iter + 1, npar))
  all.lik <- rep(NA, N.iter + 1)
  chain[1,] = start.value
  all.lik[1] <- posterior.fun(start.value, ...)
  
  for (i in 1:N.iter){
    if(i / 1000 == floor(i / 1000)){
      if(verbose) print(i)
      #chain_last1000 <- data.frame(chain[(i - 1000 + 1):i, ])
      if(i >= 10000) chain_last10000 <- chain[(i - 10000 + 1):i, ] else chain_last10000 <- chain[1:i, ]
      chain_last10000 <- data.frame(chain_last10000)
      acceptance <- 1 - mean(duplicated(chain_last10000)) # assess acceptance rate
      if(verbose) cat("acceptance rate in last 10000 generations: ", acceptance, "\n")
      if(verbose) cat("likelihood: ", all.lik[i], "\n")
      #chain_last1000_notduplicated <- chain_last1000[!duplicated(chain_last1000),]
      #n_unique_1000 <- nrow(chain_last1000_notduplicated) # number parameter sampled over last 1000
      chain_last10000_notduplicated <- chain_last10000[!duplicated(chain_last10000), ]    # unique values
      n_unique_10000 <- nrow(chain_last10000_notduplicated)
      print(n_unique_10000)
      
      if(n_unique_10000 == 1) {
        print("no new parameter accepted; reducing the variance of the proposal distribution")
        proposal.cov <- proposal.cov / factor.cov
      } else {
        if(npar > 1 & acceptance >= 0.05){
          proposal.cov <- cov(chain_last10000_notduplicated)/(10 * factor.cov)                  # update proposal cov
        }
        if(npar == 1 & acceptance >= 0.05){
          proposal.cov <- var(chain_last10000_notduplicated)/(10 * factor.cov)
        }
      }
      if(acceptance < 0.05 & i >= 10000 & n_unique_10000 != 1){   # decrease factor if acceptance rate is too small
        proposal.cov <- proposal.cov / factor.cov
        print("acceptance rate too small, reducing covariance matrix to:")
        print(diag(proposal.cov))
      }
      if(acceptance > 0.5 & i >= 10000){                         # increase factor if acceptance rate is too large
        proposal.cov <- proposal.cov * factor.cov
        print("acceptance rate too large, increasing covariance matrix to:")
        print(diag(proposal.cov))
      }
      if(verbose) print("updating proposal covariance...")
      if(verbose) cat("diagonal of variance-covariance matrix is now: ", diag(proposal.cov/factor.cov), "\n")
      print("last parameter values: ")
      print(chain[i, ])
    }
    # draw new parameters
    if(npar > 1){
      proposal <- chain[i,] + mvrnorm(1, mu = rep(0, npar), Sigma = proposal.cov) # multivariate normal distribution
    } else {
      proposal <- chain[i,] + rnorm(1, mean = 0, sd = sqrt(proposal.cov)) # univariate normal
    }
    # NOTE HERE WE CANNOT SIMPLY DRAW NEW PROPOSAL UNTIL WE FIND SUITABLE ONE - THIS WOULD BIAS THE POSTERIOR SAMPLE https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
    # see also https://stats.stackexchange.com/questions/73885/mcmc-on-a-bounded-parameter-space
    
    if(all(proposal > lower & proposal < upper)) {
      all.lik[i + 1] <- posterior.fun(proposal, ...)
    } else {
      all.lik[i + 1] <- NA # set to NA if unfeasible value
    }
    
    if(!is.na(all.lik[i+1])) probab <- exp(all.lik[i + 1] - all.lik[i]) else probab <- 0 # and set acceptance proba to 0 if NA...
    
    if(runif(1) < probab){
      chain[i+1,] <- proposal
    } else{
      chain[i+1,] <- chain[i,]
      all.lik[i + 1] <- all.lik[i]
    }
  }
  return(list(chain = chain, all.lik = all.lik))
}

nmkb <- function (par, fn, lower = -Inf, upper = Inf, control = list(), mytol = 1e-2, ...)  # from dfoptim package
{
  ctrl <- list(tol = mytol, maxfeval = min(5000, max(1500, 
                                                     20 * length(par)^2)), regsimp = TRUE, maximize = FALSE, 
               restarts.max = 3, trace = FALSE)
  namc <- match.arg(names(control), choices = names(ctrl), 
                    several.ok = TRUE)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  if (!is.null(names(control))) 
    ctrl[namc] <- control
  ftol <- ctrl$tol
  maxfeval <- ctrl$maxfeval
  regsimp <- ctrl$regsimp
  restarts.max <- ctrl$restarts.max
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  n <- length(par)
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1])/(upper[c1] - 
                                               lower[c1]) - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + 
                                                          tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  if (length(lower) == 1) 
    lower <- rep(lower, n)
  if (length(upper) == 1) 
    upper <- rep(upper, n)
  if (any(c(par <= lower, upper <= par))) 
    stop("Infeasible starting values!", call. = FALSE)
  low.finite <- is.finite(lower)
  upp.finite <- is.finite(upper)
  c1 <- low.finite & upp.finite
  c2 <- !(low.finite | upp.finite)
  c3 <- !(c1 | c2) & low.finite
  c4 <- !(c1 | c2) & upp.finite
  if (all(c2)) 
    stop("Use `nmk()' for unconstrained optimization!", call. = FALSE)
  if (maximize) 
    fnmb <- function(par) -fn(ginv(par), ...)
  else fnmb <- function(par) fn(ginv(par), ...)
  x0 <- g(par)
  if (n == 1) 
    stop(call. = FALSE, "Use `optimize' for univariate optimization")
  if (n > 30) 
    warning("Nelder-Mead should not be used for high-dimensional optimization")
  V <- cbind(rep(0, n), diag(n))
  f <- rep(0, n + 1)
  f[1] <- fnmb(x0)
  V[, 1] <- x0
  scale <- max(1, sqrt(sum(x0^2)))
  if (regsimp) {
    alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, 
                                     sqrt(n + 1) - 1)
    V[, -1] <- (x0 + alpha[2])
    diag(V[, -1]) <- x0[1:n] + alpha[1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  else {
    V[, -1] <- x0 + scale * V[, -1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  f[is.nan(f)] <- Inf
  nf <- n + 1
  ord <- order(f)
  f <- f[ord]
  V <- V[, ord]
  rho <- 1
  gamma <- 0.5
  chi <- 2
  sigma <- 0.5
  conv <- 1
  oshrink <- 1
  restarts <- 0
  orth <- 0
  dist <- f[n + 1] - f[1]
  v <- V[, -1] - V[, 1]
  delf <- f[-1] - f[1]
  diam <- sqrt(colSums(v^2))
  sgrad <- c(solve(t(v), delf))
  alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
  simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 
                                                              1])))
  itc <- 0
  conv <- 0
  message <- "Succesful convergence"
  while (nf < maxfeval & restarts < restarts.max & dist > ftol & 
         simplex.size > 1e-06) {
    fbc <- mean(f)
    happy <- 0
    itc <- itc + 1
    xbar <- rowMeans(V[, 1:n])
    xr <- (1 + rho) * xbar - rho * V[, n + 1]
    fr <- fnmb(xr)
    nf <- nf + 1
    if (is.nan(fr)) 
      fr <- Inf
    if (fr >= f[1] & fr < f[n]) {
      happy <- 1
      xnew <- xr
      fnew <- fr
    }
    else if (fr < f[1]) {
      xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 
                                                     1]
      fe <- fnmb(xe)
      if (is.nan(fe)) 
        fe <- Inf
      nf <- nf + 1
      if (fe < fr) {
        xnew <- xe
        fnew <- fe
        happy <- 1
      }
      else {
        xnew <- xr
        fnew <- fr
        happy <- 1
      }
    }
    else if (fr >= f[n] & fr < f[n + 1]) {
      xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, 
                                                       n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc <= fr) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    else if (fr >= f[n + 1]) {
      xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc < f[n + 1]) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    if (happy == 1 & oshrink == 1) {
      fbt <- mean(c(f[1:n], fnew))
      delfb <- fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) 
          cat("Trouble - restarting: \n")
        restarts <- restarts + 1
        orth <- 1
        diams <- min(diam)
        sx <- sign(0.5 * sign(sgrad))
        happy <- 0
        V[, -1] <- V[, 1]
        diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
      }
    }
    if (happy == 1) {
      V[, n + 1] <- xnew
      f[n + 1] <- fnew
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    else if (happy == 0 & restarts < restarts.max) {
      if (orth == 0) 
        orth <- 1
      V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
      for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
      nf <- nf + n
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
    f[is.nan(f)] <- Inf
    dist <- f[n + 1] - f[1]
    sgrad <- c(solve(t(v), delf))
    if (trace & !(itc%%2)) 
      cat("iter: ", itc, "\n", "value: ", f[1], "\n")
  }
  if (dist <= ftol | simplex.size <= 1e-06) {
    conv <- 0
    message <- "Successful convergence"
  }
  else if (nf >= maxfeval) {
    conv <- 1
    message <- "Maximum number of fevals exceeded"
  }
  else if (restarts >= restarts.max) {
    conv <- 2
    message <- "Stagnation in Nelder-Mead"
  }
  return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, 
              feval = nf, restarts = restarts, convergence = conv, 
              message = message))
}
myattach <- function(mylist){
  stopifnot(is.list(mylist))
  mynames <- names(mylist)
  n <- length(mylist)
  stopifnot(length(mynames) == n)
  cat("attaching: ", mynames, "\n")
  for(i in 1:n) assign(x = mynames[i], value = mylist[[i]], envir = .GlobalEnv)
}




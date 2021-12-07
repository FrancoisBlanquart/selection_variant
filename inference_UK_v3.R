rm(list = ls())

library(MASS)
library(plyr)
library(VGAM)
library(plyr)
library(Rcpp)
library(numDeriv)
library(RColorBrewer)

source("~/ownCloud/coronavirus/variant_France/clean_code/functions_v3.R")

source("~/ownCloud/coronavirus/variant_France/clean_code/get_clean_data.R")

# INFERENCE FUNCTIONS
source("~/ownCloud/coronavirus/variant_France/clean_code/inference_functions.R")

# CONCEPTUAL FIGURES
source("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figures_v2.R")

get_initRW <- function(sim_tab, varcov_mat, fixed_s3, initial_ta = NULL, initial_da = 0, weekly){
  
  # A function to generate a good starting point for RW
  # Create a good initial RW vector
  # Based on rH and rV and a priori idea of transmission advantage and duration advantage
  # And weigthed by the inverse variances
  
  if(!weekly){
    weekly_factor <- 1 
    nrows <- max_time
  } else {
    weekly_factor <- 7
    burnin <- 1
    nrows <- nrow(sim_tab) + burnin # the burnin period was already included in nweeks = nrows(sim_tab)
  }
  initRW1 <- rtoR(sim_tab$rH_error[2:(nrows-burnin)], MUW/weekly_factor, SDW/weekly_factor) # initial values for RW (using data)
  initRW2 <-  rtoR(sim_tab$rV_error[2:(nrows-burnin)], MUW * (1 + initial_da)/weekly_factor, SDW * (1 + fixed_s3)/weekly_factor)
  
  initRW1[is.na(initRW1)] <- 0.5 # too small values of r lead to NA
  initRW1[initRW1 > 5] <- 5
  initRW2[is.na(initRW2)] <- 0.5
  initRW2[initRW2 > 5] <- 5
  
  # estimate transmission advantage if no value was provided
  if(is.null(initial_ta)) initial_ta <- median(1 - initRW1/initRW2, na.rm = T)
  if(initial_ta < 0.1) initial_ta <- 0.1
  initRW2 <- initRW2 / (1+initial_ta)
  
  # weights:
  varH <- diag(varcov_mat)[seq(1, 2*(nrows-burnin)-2, 2)]; NH <- 1/varH
  varV <- diag(varcov_mat)[seq(2, 2*(nrows-burnin)-2, 2)]; NV <- 1/varV
  stopifnot(length(NH) == length(initRW1))
  stopifnot(length(NH) == length(initRW2))
  initRW <- (NH * initRW1 + NV * initRW2) / (NH + NV)
  
  return(list(initRW = initRW, initial_ta = initial_ta, initial_da = initial_da))
  
}

######################################################################################################################## 
##############################            AN ATTEMPT at A LIKELIHOOD APPROACH             ##############################
######################################################################################################################## 

# check the correlations between rV and rH
summary(lm_uk <- lm(rV ~ rH, data = uk))
summary(lm_london <- lm(rV ~ rH, data= london))
summary(lm_uk2 <- lm(rV ~ rH, data = uk2))

# get dataset:
# tab <- uk[, ]
# tab <- london[london$week_number <= 52, ]

change_names <- function(tab){
  stopifnot(all(c("N", "rH", "rV", "NVOC", "weekly_cases") %in% names(tab)))
  names(tab)[names(tab) == "N"] <- "Ntrue"
  names(tab)[names(tab) == "rH"] <- "rH_error"
  names(tab)[names(tab) == "rV"] <- "rV_error"
  names(tab)[names(tab) == "weekly_cases"] <- "Psmoothed"
  tab$fVOC <- tab$NVOC / tab$Ntrue
  stopifnot(all(c("Psmoothed", "Ntrue", "fVOC", "rH_error", "rV_error") %in% names(tab)))
  if(any(is.infinite(tab$rH_error)) | any(is.infinite(tab$rV_error))) stop()
  return(tab)
}

do_inference <- function(tab, return_full = F, varcov_weekly = NULL, fixed_s3, initial_da_list = NULL, n.iter = 5){
  
  # change column names:
  tab <- change_names(tab)
  
  # Var-covar matrix of the multivariate normal distribution
  if(is.null(varcov_weekly)) varcov_weekly <- create_cov_mat(tab = tab)
  stopifnot(isSymmetric.matrix(varcov_weekly))
  iv <- solve(varcov_weekly)
  logdetvars <- as.numeric(determinant(varcov_weekly, logarithm = T)$modulus)
  
  # data
  nweeks <- nrow(tab)
  data <- c(rbind(tab$rH_error[2:nweeks], tab$rV_error[2:nweeks])) # estimated r; what we are trying to fit; these data are in the order rH1, rV1, rH2, rV2, rH3, rV3, etc.
  #initRW <- rtoR(tab$rH_error[2:nweeks], MUW/7, SDW/7) # initial values for RW (using data)
  
  opt2_list <- list()
  all_initRW <- c()
  idx_list <- 1
  
  if(is.null(initial_da_list)) initial_da_list <- c(-0.4, -0.1, -0.05, 0, 0.05, 0.1, 0.4)
  
  for(my_initial_da in initial_da_list){
    
    # generate a good initRW:
    tmp <- get_initRW(sim_tab = tab, varcov_mat = varcov_weekly, fixed_s3, initial_ta = NULL, initial_da = my_initial_da, weekly = TRUE)
    initRW <- tmp$initRW; all_initRW <- rbind(all_initRW, tmp$initRW)
    initial_ta <- tmp$initial_ta
    initial_da <- tmp$initial_da
    if(any(is.na(initRW))) stop("Some elements of initRW are NA")
    
    lower2 <- c(0, -0.8)
    upper2 <- c(3, 3)
    npar2 <- 2
    
    initfun2 <- function() c(initial_ta, initial_da)
    opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = fastneglik2, init.fun = initfun2, optim.fun = "optim", lower = lower2, upper = upper2, mydata = data, myinitRW = initRW, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = fixed_s3, weekly = TRUE)
    opt2_list[[idx_list]] <- opt2; idx_list <- idx_list + 1
    #2nd round:
    
    for(l in 1:n.iter){
      tmp <- get_initRW(sim_tab = tab, varcov_mat = varcov_weekly, fixed_s3, initial_ta = opt2$pars[1], initial_da = opt2$pars[2], weekly = TRUE)
      initRW <- (tmp$initRW); all_initRW <- rbind(all_initRW, tmp$initRW)
      initial_ta <- tmp$initial_ta
      initial_da <- tmp$initial_da
      if(initial_da >= upper2[2]) initial_da <- upper2[2] * 0.9
      if(initial_ta >= upper2[1]) initial_ta <- upper2[1] * 0.9
      initfun2 <- function() c(initial_ta, initial_da)
      opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = fastneglik2, init.fun = initfun2, optim.fun = "optim", lower = lower2, upper = upper2, mydata = data, myinitRW = initRW, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = fixed_s3, weekly = TRUE)
      opt2_list[[idx_list]] <- opt2; idx_list <- idx_list+1
    }
  }
  
  all_lik <- unlist(lapply(opt2_list, function(x) x$lik))
  stopifnot(length(opt2_list) == length(all_lik))
  stopifnot(length(opt2_list) == nrow(all_initRW))
  opt2 <- opt2_list[[which.min(all_lik)]]
  initRW <- all_initRW[which.min(all_lik), ]
  
  #plot(all_lik, type = "o", pch = 20)
  #plot(initRW, type = "l", lwd = 3, pch= 20, col= "blue", ylim = c(0,2))
  
  if(return_full){
    return(list(
      opt = opt2, mydata = data, myinitRW = initRW, vars = varcov_weekly
      ))
  } else {
    return(list(opt = opt2))
  }
}

# subsets of data used for inference:
sub_uk <- uk$week_number <= 62
sub_london <- london$week_number <= 62
sub_uk2 <- uk2$day <= 530

########################################################################### INFERENCE FOR ALPHA VARIANT ###########################################################################

# 1. Heuristic optimisation
# 1.1 Exploring different reduction in weekly cases to account for negative binomial response

liks <- c()
factors <- 1.1^(0:20)
for(ff in factors){
  uk_eff <- uk[sub_uk, ]
  uk_eff$weekly_cases <- round(uk_eff$weekly_cases/ff)
  opt_uk <- do_inference(tab =  uk_eff, return_full = T, fixed_s3 = 0, n.iter = 20)
  print(opt_uk$opt)
  liks <- c(liks, opt_uk$opt$lik)
}
plot(factors, liks, type = "l", las = 1, bty= "n", xlab = c("Factor reducing effective number of cases", ylab = "Likelihood"))
best_factor <- factors[which.min(liks)] # 2.45 in previous iteration


# 1.2 redo inference with best effective number of cases:
uk_eff <- uk[sub_uk, ]
uk_eff$weekly_cases <- (uk_eff$weekly_cases/best_factor)
opt_uk <- do_inference(tab =  uk_eff, return_full = T, fixed_s3 = 0, n.iter = 100)

opt_uk$opt
#ci_uk$boot_ci

# 2. now full optimisation (inferring all the R0)

lower0 <- c(0.3, -0.5, rep(0.4, sum(sub_uk) - 1))
upper0 <- c(0.7, 1, rep(1.6, sum(sub_uk) - 1))
init_fun0 <- function() c(runif(n = 1, min = lower0[1], max = upper0[1]), runif(1, min = lower0[2], max = upper0[2]), opt_uk$myinitRW)

# 2.1. tune the error variance again with full optimisation
# the choice of factor intervals is based on several explorations
liks <- c()
factors1 <- 10^seq(0, 1, 0.1)
factors2 <- 10^seq(2, 4, 0.2)
all_opt_uk_full <- list()
idx <- 1
for(f1 in factors1){
  for(f2 in factors2){
    print(c(f1, f2))
    uk_eff <- uk[sub_uk, ]
    uk_eff$weekly_cases <- uk_eff$weekly_cases/f1 # reduced effective number of cases
    uk_eff$NVOC <- uk_eff$NVOC / f2; uk_eff$N <- uk_eff$N / f2 # reduced effective cases assayed for VOC
    varcov_weekly <- create_cov_mat(tab = change_names(uk_eff))
    iv <- solve(varcov_weekly)
    logdetvars <- as.numeric(determinant(varcov_weekly, logarithm = T)$modulus)
    nweeks <- nrow(uk_eff)
    mydata <- c(t(uk_eff[2:nweeks, c("rH", "rV")]))
    opt_uk_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                      mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)

    print(opt_uk_full)
    all_opt_uk_full[[idx]] <- opt_uk_full; idx <- idx + 1
    liks <- rbind(liks, c(f1, f2, opt_uk_full$lik))
  }
}
plot(NULL, xlim = c(0,10000), ylim = c(-100, -10), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Log likelihood")
for(ii in 1:11){
  sub <- which(liks[, 1] == factors1[ii])
  points(liks[sub, 2], liks[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
which.min(liks[, 3]) -> best_idx
best_f1 <- liks[best_idx, 1]
best_f2 <- liks[best_idx, 2]

# fastneglik(allpars = init_fun0(), mydata = mydata, inverse_vars = solve(varcov_weekly), logdetvars = as.numeric(determinant(varcov_weekly, logarithm = T)$modulus), fixed_s3 = 0, weekly = T)

# 2.2 full optimisation with the best error variance

uk_eff <- uk[sub_uk, ]
uk_eff$weekly_cases <- uk_eff$weekly_cases / best_f1 # reduced effective number of cases
uk_eff$NVOC <- uk_eff$NVOC / best_f2; uk_eff$N <- uk_eff$N / best_f2 # reduced effective cases assayed for VOC

varcov_weekly <- create_cov_mat(tab = change_names(uk_eff))
iv <- solve(varcov_weekly)
logdetvars <- as.numeric(determinant(varcov_weekly, logarithm = T)$modulus)
nweeks <- nrow(uk_eff)
mydata <- c(t(uk_eff[2:nweeks, c("rH", "rV")]))
init_fun0 <- function() c(runif(n = 1, min = lower0[1], max = upper0[1]), runif(1, min = lower0[2], max = upper0[2]), opt_uk_full$pars[3:(nweeks+1)])
opt_uk_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)

h <- hessian(fastneglik, x = opt_uk_full$pars,  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T) # Hessian at the optimum
cov_matrix <- solve(h) # covariance matrix is inverse of Hessian https://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s?rq=1
sds <- sqrt(diag(cov_matrix))
ci_hessian <- rbind(opt_uk_full$pars - 1.96 * sds, opt_uk_full$pars, opt_uk_full$pars + 1.96 * sds)
rownames(ci_hessian) <- c("2.5%", "ML", "97.5%")

plot(NULL, xlim = c(0, 30), ylim = c(0.4, 1.6), bty = "n", ylab = "R", xlab = "Week", las = 1)
#points(opt_uk$myinitRW, type = "o", pch = 20, col = "red") # Heuristic is not so bad
polygon(
  c(1:28, rev(1:28)),
  c(ci_hessian["2.5%",3:30], rev(ci_hessian["97.5%",3:30])),
  col = adjustcolor("red", 0.5), border = NA
)
points(opt_uk_full$pars[3:(sum(sub_uk)+2)], type = "o", pch = 20)


########################################################################### TEMPORAL INFERENCE FOR DELTA VARIANT ###########################################################################


varcov_weekly <- create_cov_mat(tab = change_names(uk2[sub_uk2, ]))
myinitRW <- get_initRW(sim_tab = change_names(uk2[sub_uk2, ]), varcov_mat = varcov_weekly, fixed_s3 = 0, weekly = TRUE)$initRW

# 2. now full optimisation (inferring all the R0)

lower0 <- c(0.3, -0.5, rep(0.4, sum(sub_uk2) - 1))
upper0 <- c(4, 1, rep(1.6, sum(sub_uk2) - 1))
init_fun0 <- function() c(runif(n = 1, min = lower0[1], max = upper0[1]), runif(1, min = lower0[2], max = upper0[2]), myinitRW)

# 2.1. tune the error variance again with full optimisation
# the choice of factor intervals is based on several explorations
liks2 <- c()
factors22 <- 10^seq(0, 2, 0.5)
all_opt_uk2_full <- list()
idx <- 1
for(f1 in factors1){
  for(f2 in factors22){
    print(c(f1, f2))
    uk2_eff <- uk2[sub_uk2, ]
    uk2_eff$weekly_cases <- uk2_eff$weekly_cases/f1 # reduced effective number of cases
    uk2_eff$NVOC <- uk2_eff$NVOC / f2; uk2_eff$N <- uk2_eff$N / f2 # reduced effective cases assayed for VOC
    varcov_weekly <- create_cov_mat(tab = change_names(uk2_eff))
    iv <- solve(varcov_weekly)
    logdetvars <- as.numeric(determinant(varcov_weekly, logarithm = T)$modulus)
    nweeks <- nrow(uk2_eff)
    mydata <- c(t(uk2_eff[2:nweeks, c("rH", "rV")]))
    opt_uk2_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                      mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)
    
    print(opt_uk2_full)
    all_opt_uk2_full[[idx]] <- opt_uk2_full; idx <- idx + 1
    liks2 <- rbind(liks2, c(f1, f2, opt_uk2_full$lik))
  }
}
plot(NULL, xlim = c(0,100), ylim = c(-50, 200), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Log likelihood")
for(ii in 1:11){
  sub <- which(liks2[, 1] == factors1[ii])
  points(liks2[sub, 2], liks2[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
which.min(liks2[, 3]) -> best_idx
best2_f1 <- liks2[best_idx, 1]
best2_f2 <- liks2[best_idx, 2]

# 2.2 full optimisation with the best error variance

uk2_eff <- uk2[sub_uk2, ]
uk2_eff$weekly_cases <- uk2_eff$weekly_cases / best2_f1 # reduced effective number of cases
uk2_eff$NVOC <- uk2_eff$NVOC / best2_f2; uk2_eff$N <- uk2_eff$N / best2_f2 # reduced effective cases assayed for VOC

varcov_weekly <- create_cov_mat(tab = change_names(uk2_eff))
iv <- solve(varcov_weekly)
logdetvars <- as.numeric(determinant(varcov_weekly, logarithm = T)$modulus)
nweeks <- nrow(uk2_eff)
mydata <- c(t(uk2_eff[2:nweeks, c("rH", "rV")]))
init_fun0 <- function() c(runif(n = 1, min = lower0[1], max = upper0[1]), runif(1, min = lower0[2], max = upper0[2]), opt_uk2_full$pars[3:(nweeks+1)])
opt_uk2_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)

h2 <- hessian(fastneglik, x = opt_uk2_full$pars,  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T) # Hessian at the optimum
cov_matrix2 <- solve(h2) # covariance matrix is inverse of Hessian https://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s?rq=1
sds2 <- sqrt(diag(cov_matrix2))
ci2_hessian <- rbind(opt_uk2_full$pars - 1.96 * sds2, opt_uk2_full$pars, opt_uk2_full$pars + 1.96 * sds2)
rownames(ci2_hessian) <- c("2.5%", "ML", "97.5%")

################################################################################### SPATIAL INFERENCE FOR DELTA VARIANT ####################################################################################

rm(mydata)
eu <- read.csv(file = "~/ownCloud/coronavirus/variant_France/clean_code/spatial_variation_EU.csv")
# need to add a fake country so we can use same functions as for temporal data (where there is no growth rate in the first row)
eu <- rbind(eu, c("aaa", rep(NA, 9)))
eu <- eu[order(eu$country), ]
for(cc in names(eu)[2:9]) eu[, cc] <- as.numeric(as.character(eu[, cc]))

eu$fVOC0 <- eu$NVOC0/eu$Ntrue0
eu$fVOC1 <- eu$NVOC1/eu$Ntrue1
varcov_spatial <- create_spatial_cov_mat(tab = eu)

# initial RW:
names(eu)[names(eu) == "rH"] <-"rH_error"
names(eu)[names(eu) == "rV"] <-"rV_error"
myinitRW <- get_initRW(sim_tab = eu, varcov_mat = varcov_spatial, fixed_s3 = 0, weekly = TRUE)$initRW
nweeks <- nrow(eu)
mydata <- c(t(eu[2:nweeks, c("rH_error", "rV_error")]))

# 2. now full optimisation (inferring all the R0)

# 2.1. tune the error variance again with full optimisation
# the choice of factor intervals is based on several explorations
liks3 <- c()
all_opt_eu_full <- list()
idx <- 1
factors23 <- 10^seq(0, 1, 0.25)
for(f1 in factors1){
  for(f2 in factors23){
    print(c(f1, f2))
    eu_eff <- eu
    eu_eff$Psmoothed0 <- eu_eff$Psmoothed0/f1; eu_eff$Psmoothed1 <- eu_eff$Psmoothed1/f1 # reduced effective number of cases
    eu_eff$NVOC0 <- eu_eff$NVOC0 / f2; eu_eff$Ntrue0 <- eu_eff$Ntrue0 / f2; eu_eff$NVOC1 <- eu_eff$NVOC1 / f2; eu_eff$Ntrue1 <- eu_eff$Ntrue1 / f2 # reduced effective cases assayed for VOC
    varcov_spatial <- create_spatial_cov_mat(tab = eu_eff)
    iv <- solve(varcov_spatial)
    logdetvars <- as.numeric(determinant(varcov_spatial, logarithm = T)$modulus)
    nweeks <- nrow(eu_eff)
    opt_eu_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                      mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)
    
    print(opt_eu_full)
    all_opt_eu_full[[idx]] <- opt_eu_full; idx <- idx + 1
    liks3 <- rbind(liks3, c(f1, f2, opt_eu_full$lik))
  }
}
plot(NULL, xlim = c(0,10), ylim = c(-40, 10), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Log likelihood")
for(ii in 1:11){
  sub <- which(liks3[, 1] == factors1[ii])
  points(liks3[sub, 2], liks3[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
which.min(liks3[, 3]) -> best_idx
best3_f1 <- liks3[best_idx, 1]
best3_f2 <- liks3[best_idx, 2]

# 2.2 and reoptimise with best error variances
eu_eff <- eu
eu_eff$Psmoothed0 <- eu_eff$Psmoothed0 / best3_f1; eu_eff$Psmoothed1 <- eu_eff$Psmoothed1 / best3_f1 # reduced effective number of cases
eu_eff$NVOC0 <- eu_eff$NVOC0 / best3_f2; eu_eff$Ntrue0 <- eu_eff$Ntrue0 / best3_f2; eu_eff$NVOC1 <- eu_eff$NVOC1 / best3_f2; eu_eff$Ntrue1 <- eu_eff$Ntrue1 / best3_f2 # reduced effective cases assayed for VOC
varcov_spatial <- create_spatial_cov_mat(tab = eu_eff)
iv <- solve(varcov_spatial)
logdetvars <- as.numeric(determinant(varcov_spatial, logarithm = T)$modulus)
opt_eu_full <- optim.fun.repeated(n.repeats = 30, lik.fun = fastneglik, init.fun = init_fun0, optim.fun = "optim", lower = lower0, upper = upper0,
                                  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T)
opt_eu_full

# confidence intervals
h3 <- hessian(fastneglik, x = opt_eu_full$pars,  mydata = mydata, inverse_vars = iv, logdetvars = logdetvars, fixed_s3 = 0, weekly = T) # Hessian at the optimum
cov_matrix3 <- solve(h3) # covariance matrix is inverse of Hessian https://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s?rq=1
sds3 <- sqrt(diag(cov_matrix3))
ci3_hessian <- rbind(opt_eu_full$pars - 1.96 * sds3, opt_eu_full$pars, opt_eu_full$pars + 1.96 * sds3)
rownames(ci3_hessian) <- c("2.5%", "ML", "97.5%")


save.image("~/ownCloud/coronavirus/variant_France/clean_code/inference_results_v2.RData")


################################################################################### PLOT ####################################################################################
rm(list =ls())
load("~/ownCloud/coronavirus/variant_France/clean_code/inference_results_v2.RData")
library(RColorBrewer)

pdf("~/ownCloud/coronavirus/variant_France/clean_code/inferred_R_t.pdf", width = 6, height = 4)
par(mfrow = c(1,1), mar = c(4,4,1,1))

plot(NULL, xlim = c(32, 75), ylim = c(0.4, 1.6), bty = "n", ylab = "R(t)", xlab = "Week number", las = 1, axes = F)
axis(1, seq(35, 75, 5))
axis(2, seq(0.4, 1.6, 0.2), las = 1)

#points(opt_uk$myinitRW, type = "o", pch = 20, col = "red") # Heuristic is not so bad
polygon(
  c(35:62, 62:35),
  c(ci_hessian["2.5%",3:30], rev(ci_hessian["97.5%",3:30])),
  col = adjustcolor("red", 0.5), border = NA
)
points(35:62, ci_hessian["ML",3:30], type = "o", pch = 20)

polygon(
  c(64:75, 75:64),
  c(ci2_hessian["2.5%",3:14], rev(ci2_hessian["97.5%",3:14])),
  col = adjustcolor("red", 0.5), border = NA
)
points(64:75, ci2_hessian["ML",3:14], type = "o", pch = 20)

text(x = 49, y = 1.5, "(1) Alpha variant")
text(x = 70, y = 1.5, "(2) Delta variant", cex = 1)
abline(v = 63, lty = 2)
dev.off()


pdf("~/ownCloud/coronavirus/variant_France/clean_code/tuning_error_variance.pdf", width = 9, height = 6)
par(mfrow = c(2,3), mar = c(4,4,1,1))
# likelihood 
plot(NULL, xlim = c(0,10000), ylim = c(-100, -10), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Negative log likelihood")
for(ii in 1:11){
  sub <- which(liks[, 1] == factors1[ii])
  points(liks[sub, 2], liks[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks[which.min(liks[, 3]), 2], liks[which.min(liks[, 3]), 3], cex = 1.5)

plot(NULL, xlim = c(0,100), ylim = c(-50, 0), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Negative log likelihood")
for(ii in 1:11){
  sub <- which(liks2[, 1] == factors1[ii])
  points(liks2[sub, 2], liks2[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks2[which.min(liks2[, 3]), 2], liks2[which.min(liks2[, 3]), 3], cex = 1.5)

plot(NULL, xlim = c(0,10), ylim = c(-30, -5), las = 1, bty = "n", xlab = "Factor reducing cases assayed for VOC", ylab = "Negative log likelihood")
for(ii in 1:11){
  sub <- which(liks3[, 1] == factors1[ii])
  points(liks3[sub, 2], liks3[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks3[which.min(liks3[, 3]), 2], liks3[which.min(liks3[, 3]), 3], cex = 1.5)

plot(NULL, xlim = c(0,10), ylim = c(-100, -10), las = 1, bty = "n", xlab = "Factor reducing cases number", ylab = "Megative log likelihood")
for(ii in 1:11){
  sub <- which(liks[, 2] == factors2[ii])
  points(liks[sub, 1], liks[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks[which.min(liks[, 3]), 1], liks[which.min(liks[, 3]), 3], cex = 1.5)

plot(NULL, xlim = c(0,10), ylim = c(-50, 100), las = 1, bty = "n", xlab = "Factor reducing cases number", ylab = "Negative log likelihood")
for(ii in 1:5){
  sub <- which(liks2[, 2] == factors22[ii])
  points(liks2[sub, 1], liks2[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks2[which.min(liks2[, 3]), 1], liks2[which.min(liks2[, 3]), 3], cex = 1.5)

plot(NULL, xlim = c(0,10), ylim = c(-50, 0), las = 1, bty = "n", xlab = "Factor reducing cases number", ylab = "Negative log likelihood")
for(ii in 1:5){
  sub <- which(liks3[, 2] == factors23[ii])
  points(liks3[sub, 1], liks3[sub, 3], pch = 20, type = "o", col = brewer.pal(n = 12, name = "Paired")[ii])
}
points(liks3[which.min(liks3[, 3]), 1], liks3[which.min(liks3[, 3]), 3], cex = 1.5)

dev.off()

# plot rH-rV relations

plot_ip <- function(s1, s2, s3, sdvec, mycol, ...){ # ML parameters and bootstrapped parameters
  
  nbs <- 1000
  x <- seq(0, 15, 0.01)
  boot_s1 <- s1 + rnorm(nbs, mean = 0, sd = sdvec[1])
  boot_s2 <- s2 + rnorm(nbs, mean = 0, sd = sdvec[2])
  
  gamma_pars_W <- get_gamma_function(mu = MUW            , sd = SDW)
  gamma_pars_M <- get_gamma_function(mu = MUW * (1 + s2) , sd = SDW * (1 + s3))
  
  # bootstrapped parameters:
  gamma_pars_W_boot <- get_gamma_function(mu = MUW            , sd = SDW)
  gamma_pars_M_boot <- get_gamma_function(mu = MUW * (1 + boot_s2) , sd = SDW * (1 + s3))
  
  # plot:
  plot(NULL, xlim = c(0, 15), ylim = c(0, 0.5), bty = "n", ylab = "Infectiousness profile", xlab = "Days", las = 1, ...)
  
  tmp <- sapply(1:nbs, function(i) (1 + boot_s1[i]) * dgamma(x, shape = gamma_pars_M_boot$shape[i], scale = gamma_pars_M_boot$scale[i]))
  tmp <- apply(tmp, 1, quantile, c(0.025, 0.975), na.rm = T)
  polygon(c(x, rev(x)), y = c(tmp["2.5%", ], rev(tmp["97.5%", ])), col = adjustcolor(mycol, 0.5), border = NA)
  
  points(x <- seq(0, 15, 0.1), dgamma(x, shape = gamma_pars_W$shape, scale = gamma_pars_W$scale), type = "l", lwd = 3)
  points(x, (1 + s1) * dgamma(x, shape = gamma_pars_M$shape, scale = gamma_pars_M$scale), type = "l", col = mycol, lwd = 3)
  
}
plot_contour <- function(pars, covmat, mycol){
   # test
  # pars = opt_uk_full$pars; covmat = cov_matrix; mycol = cols[1]
  
  n <- length(pars)
  stopifnot(all(dim(covmat) == c(n, n)))
  boot_pars <- MASS::mvrnorm(n = 1000, mu = pars, Sigma = covmat)
  boot_means <- apply(boot_pars, 1, function(x) getmeans(pars = c(x[1:2], 0), Rpars = x[3:n], weekly = T))
  boot_meansH <- boot_means[seq(1, (n-2)*2, 2), ]; boot_meansV <- boot_means[seq(2, (n-2)*2, 2), ]
  
  means <- getmeans(pars = c(pars[1:2], 0), Rpars = pars[3:n], weekly = T)
  meansH <- means[seq(1, (n-2)*2, 2)]; meansV <- means[seq(2, (n-2)*2, 2)]
  
  # border where we want to infer CI:
  minx <- min(meansH)
  maxx <- max(meansH)
  seq_x <- seq(minx, maxx, (maxx-minx)/100) # sequences of x for which we want a CI
  
  all_y <- matrix(NA, ncol = 101, nrow = 1000)
  for(i in 1:1000){
    xH <- boot_meansH[, i]
    xV <- boot_meansV[, i]
    if(all(!is.na(xV))){
      coeff <- lm(xV ~ xH)$coefficients # fit a linear model to the xH, xV relationship to smooth it and have smooth CI
      all_y[i, ] <- seq_x *coeff["xH"]  + coeff["(Intercept)"]
    }
  }
  all_y <- apply(all_y, 2, quantile, c(0.025, 0.975), na.rm = T)
  polygon(x = c(seq_x, rev(seq_x)), y = c(all_y[1, ], rev(all_y[2, ])), col = adjustcolor(mycol, 0.5), border = NA)
  points(meansH[order(meansH)], meansV[order(meansV)], col = mycol, lwd = 3, type = "l")
  
}

pdf("~/ownCloud/coronavirus/variant_France/clean_code/figure_profiles.pdf", width = 6, height = 6)

par(mar = c(4,4,1,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE), width = c(1.5, 1), height = c(1, 1))

#xlab <- expression(paste("Effective reproduction number ", R, " wild type"))
plot(NULL, xlim = c(-0.5, 0.5), ylim  = c(-0.5, 0.5), xlab = "Growth rate H (per week)", ylab = "Growth rate E (per week)", las = 1, bty = "n")
abline(h = 0, lty = 2, col = "gray")
abline(v = 0, lty = 2, col = "gray")
abline(0,1, col = "gray")


x1 <- 7*sapply(Rfactor, function(f) with(get_mutant_par(sc1), Rtor(R = f*RW, mu = MUW, SD = SDW)))
x2 <- 7*sapply(Rfactor, function(f) with(get_mutant_par(sc2), Rtor(R = f*RW, mu = MUW, SD = SDW)))
y1 <- 7*sapply(Rfactor, function(f) with(get_mutant_par(sc1), Rtor(R = f*RV, mu = MUV, SD = SDV)))
y2 <- 7*sapply(Rfactor, function(f) with(get_mutant_par(sc2), Rtor(R = f*RV, mu = MUV, SD = SDV)))

sub1 <- which(x1 > -0.5 & x1 < 0.5 & y1 > -0.5 & y1 < 0.5)
sub2 <- which(x2 > -0.5 & x2 < 0.5 & y2 > -0.5 & y2 < 0.5)

cols2 <- c(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:4], mygray) # color scheme for conceptual sim
points(x1[sub1], y1[sub1], col = cols2[1], lty = 1, type = "l", lwd = 3)
points(x2[sub2], y2[sub2], col = cols2[2], lty = 1, type = "l", lwd = 3)

points(x = 7 * Rtor(R = 1/(1 + sc1$s1), mu = MUW, SD = SDW), y = 0, pch = 20, cex = 2)
#text(x = rep(1.85, 6), y = c(0.003, 0.021, 0.04, 0.057, 0.073, 0.09), labels = c("WT", "+10%", "+20%", "+30%", "+40%", "+50%"), col = mygray)
legend("bottomright", legend = c("Shorter mgt", "Longer mgt"), lwd = 3, col = c(cols2[1:2]), bty = "n")

segments(x0 = x2[100], y0 = y2[100], x1 = x2[120], y1 = y2[100], lwd = 1.5, col = cols2[2])
segments(x0 = x2[120], y0 = y2[100], x1 = x2[120], y1 = y2[120], lwd = 1.5, col = cols2[2])

text(x = -0.2, y = 0.3, labels = expression(paste("Slope", phantom() %~~% phantom(), "1 + ", s[1], cv[H]^2, " - ", s[2])), col = cols2[2])


cols <- RColorBrewer::brewer.pal(n = 8, name = "Paired")[c(2, 1, 8)] # colors
plot(NULL, xlim = c(-1, 1.1), ylim  = c(-1, 1.51), xlab = "Growth rate H (per week)", ylab = "Growth rate E (per week)", las = 1, bty = "n")
abline(h = 0, lty = 2, col = "gray")
abline(v = 0, lty = 2, col = "gray")
abline(0,1, col = "gray")
plot_contour(pars = opt_uk_full$pars, covmat = cov_matrix, mycol = cols[1])
#plot_contour(pars = opt_uk2_full$pars, covmat = cov_matrix2, mycol = cols[3])
plot_contour(pars = opt_eu_full$pars, covmat = cov_matrix3, mycol = cols[3])

points(uk$rH[sub_uk], uk$rV[sub_uk], pch = 20, col = cols[1], type = "p")
#points(uk2$rH, uk2$rV, pch = 20, col = cols[3], type = "p")
points(eu$rH, eu$rV, pch = 20, col = cols[3], type = "p")

legend("bottomright", legend = c("Alpha", "Delta"), col = c(cols[1], cols[3]), pch = 20, lwd =1, bty = "n", cex = 1)

# now infectiousness profiles:
plot_ip(s1 = opt_uk_full$pars[1], s2 = opt_uk_full$pars[2], s3 = 0, sdvec = sds, mycol = cols[1]) #, main  ="Alpha variant")
#plot_ip(s1 = opt_uk2_full$pars[1], s2 = opt_uk2_full$pars[2], s3 = 0, sdvec = sds2, mycol = cols[3], main = "Delta variant")
plot_ip(s1 = opt_eu_full$pars[1] * (1 + opt_uk_full$pars[1]), s2 = opt_eu_full$pars[2], s3 = 0, sdvec = sds3, mycol = cols[3]) #, main = "Delta variant")
# note the s1 of the Delta is relative to the Alpha

dev.off()


# figure showing frequencies and cases
myfactor <- 0.8
pdf("~/ownCloud/coronavirus/variant_France/clean_code/figure_inference_v2.pdf", width = 6 * 2 * myfactor, height = 4 * 2 * myfactor)

par(mar = c(6,4,1,6))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), width = c(1, 1), height = c(1, 1))

#par(mfrow = c(1,3))

plot(uk$week_number, uk$NVOC/uk$N, xlab = "Week number", ylab = "Frequency variant", las = 1, xaxs = "i", yaxs = "i", las = 1, lwd = 3, col = cols[1], type = "o", pch = 20, ylim = c(0, 1), xlim = c(34, 60), axes = F)
points(uk$week_number, uk$weekly_cases/400000, col = "gray", lwd = 3, pch = 20, type = "o")
points(uk$week_number, uk$N/400000, col = "gray", lwd = 3, pch = 20, type = "o", lty = 3)
axis(side = 1, at = seq(35, 60, 5))
axis(side = 2, at = seq(0, 1, 0.2), las = 1)
axis(side = 4, at = seq(0, 1, 0.25), labels = seq(0, 400000, 100000), las = 1)
mtext("Weekly cases number", side=4, line=4, cex = 0.7)
legend("topleft", lwd = 3, pch = 20, col = c(cols[1], "gray", "gray"), lty = c(1,1,3), legend = c("Alpha frequency", "Cases number", "Sample size"), bty = "n", cex = 0.9)

if(plot_london <- F){
  plot(london$week_number, london$NVOC/london$N, xlab = "Week number", ylab = "Frequency Alpha variant", las = 1, xaxs = "i", yaxs = "i", las = 1, lwd = 3, col = cols[2], type = "o", pch = 20, ylim = c(0, 1), xlim = c(34, 60), axes = F)
  points(london$week_number, london$weekly_cases/400000, col = "gray", lwd = 3, pch = 20, type = "o")
  points(london$week_number, london$N/400000, col = "gray", lwd = 3, pch = 20, type = "o", lty = 3)
  axis(side = 1, at = seq(35, 60, 5))
  axis(side = 2, at = seq(0, 1, 0.1), las = 1)
  axis(side = 4, at = seq(0, 1, 0.25), labels = seq(0, 400000, 100000), las = 1)
  mtext("Weekly cases number", side=4, line=4, cex = 0.7)
  legend("topleft", lwd = 3, pch = 20, col = c(cols[2], "gray", "gray"), lty = c(1,1,3), legend = c("Alpha frequency London", "Cases number", "Sample size"), bty = "n")
}

plot(uk2$week_number, uk2$NVOC/uk2$N, xlab = "Week number", ylab = "Frequency Delta variant", las = 1, xaxs = "i", yaxs = "i", las = 1, lwd = 3, col = cols[3], type = "o", pch = 20, ylim = c(0, 1), xlim = c(60, 80), axes = F)
points(uk2$week_number, uk2$weekly_cases/80000, col = "gray", lwd = 3, pch = 20, type = "o")
points(uk2$week_number, uk2$N/80000, col = "gray", lwd = 3, pch = 20, type = "o", lty = 3)
axis(side = 1, at = seq(60, 80, 5))
axis(side = 2, at = seq(0, 1, 0.2), las = 1)
axis(side = 4, at = seq(0, 1, 0.25), labels = seq(0, 80000, 20000), las = 1)
mtext("Weekly cases number", side=4, line=4, cex = 0.7)
legend("topleft", lwd = 3, pch = 20, col = c(cols[3], "gray", "gray"), lty = c(1,1,3), legend = c("Delta frequency", "Cases number", "Sample size"), bty = "n", cex = 0.9)

par(mar = c(4,23*myfactor,1,23*myfactor))
all_freq_eu <- read.csv("~/ownCloud/coronavirus/variant_France/clean_code/delta_frequencies_EU.csv")
list_cc <- c("Austria", "Belgium", "Denmark", "France", "Germany", "Greece", "Ireland", "Italy", "Netherlands", "Norway", "Sweden")
cols3 <- RColorBrewer::brewer.pal(n = 11, name = "Paired")
names(cols3) <- list_cc
plot(NULL, xlab = "Week number", ylab = "Frequency Delta variant", las = 1, xaxs = "i", yaxs = "i", las = 1, ylim = c(0, 1), xlim = c(70, 85), axes = F)
for(cc in list_cc){
  idx <- which(all_freq_eu$country==cc)
  points(all_freq_eu$week[idx] + 53, all_freq_eu$percent_variant[idx]/100, col = cols3[cc], type = "o", pch = 20, main = cc)
  
  idx_first50 <- idx[which(all_freq_eu$percent_variant[idx]/100 > 0.5)[1]]; sub <- (idx_first50-1):idx_first50
  points(all_freq_eu$week[sub] + 53, all_freq_eu$percent_variant[sub]/100, col = cols3[cc], type = "o", pch = 20, cex = 2, lwd = 3)
  
  #print(range(all_freq_eu$week[idx] + 53))
}
axis(side = 1, at = seq(60, 85, 5))
axis(side = 2, at = seq(0, 1, 0.2), las = 1)

legend("topleft", legend = names(cols3), col = cols3, pch = 20, lty = 1, bty = "n")
dev.off()

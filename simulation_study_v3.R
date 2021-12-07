rm(list = ls())

library(MASS)
library(plyr)
library(VGAM)
library(plyr)
library(Rcpp)

# SIMULATION EXPERIENCE

source("~/ownCloud/coronavirus/variant_France/clean_code/functions_v3.R")

source("~/ownCloud/coronavirus/variant_France/clean_code/get_clean_data.R")

# INFERENCE FUNCTIONS
source("~/ownCloud/coronavirus/variant_France/clean_code/inference_functions.R")

source("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figures_v2.R")

# need to load cases_per_French_region_v5.R

sourceCpp("~/ownCloud/coronavirus/variant_France/clean_code/simulate_v2.cpp")  # simulations in Rcpp

make_smooth <- function(vec){
  x <- 1:length(vec)
  y <- loess(vec ~ x, span = 0.2)
  yfit <- predict(y, newdata=  x)
  yfit
}

pop_idf <- 12278210 # pop IDF
immunity_idf <- immunity[which(immunity$date=="2021-01-01" & immunity$region == "IDF"), "infected.mean"]/100

# reset sd_advantage
sd_advantage <- 0

################################################## RUN DYNAMIC SIMULATIONS AND PLOT (RH, RV) ##################################################

# BM MOTION OF THE WT with mean 1.3 and sd 

R0H <- 1.5
flag <- T
cond_on_R0 <- function(myR0vec) mean(myR0vec[1:20]) > 1.1 & mean(myR0vec[(max_time-20):max_time]) < 0.9 & all(myR0vec > 0.5) & all(myR0vec < 5)
while(flag){
  R0vecH <- mvrnorm(n = 1, mu = rep(R0H, max_time), Sigma = create_BMautocov_mat(temporal_autocor, max_time = max_time))
  if(cond_on_R0(R0vecH)) flag <- F # select a realisation of BM where it is below 1 at the end
}
R0vecH <- c(2.1307553, 2.0105587, 2.0023232, 1.9215775, 1.9642554, 1.9238395, 2.1029669, 2.0780443, 2.1581173, 2.0758108, 2.1253340, 1.9948015, 1.8292173, 1.7309693, 1.7447241, 1.7712299, 1.7653672, 1.7129285, 1.8179133, 1.7841625, 1.9068882, 1.9366379, 2.0484530, 2.1590289, 2.1527725, 2.2736368, 2.2387613, 2.2257805, 2.1303752, 2.1895513, 2.2356056, 2.2355079, 2.3258214, 2.2258386, 2.1631581, 2.2186670, 2.1714491, 2.2069247, 2.1436307, 2.1366790, 2.0831706, 1.9912173, 1.8703819, 1.6309556, 1.4634919, 1.3115412, 1.1725386, 1.0649104, 1.0250046, 0.9133384, 1.0018016, 0.9228869, 0.9090948, 0.7077709, 0.5578058, 0.6371114, 0.7721366, 0.7062939, 0.6210061, 0.4607542, 0.5308621, 0.3945342, 0.4594047, 0.6524794, 0.6660125, 0.5297658, 0.5369410, 0.5645049, 0.6640699, 0.6662376, 0.5892468, 0.5752508, 0.7023259, 0.6168691, 0.4811847, 0.3900102, 0.4270968, 0.2734790, 0.3581384, 0.3733780)
# as of 30/08 move to smooth linear function:
Rinit <- 2; Rfin <- 0.5
R0vecH <- c(rep(Rinit, 10), seq(Rinit, Rfin, by = -(Rinit-Rfin)/(70-1)))
plot(R0vecH, type = "l"); abline(h = 1, lty = 2)

# PARAMETERS OF THE TRANSMISSION FOR VARIANTS:
burnin <- 10


beta_advantage <- 0.
sc1 <- list(s1 = beta_advantage, s2 = -0.4, s3 = 0)
sc2 <- list(s1 = beta_advantage, s2 = 0.4, s3 = 0.)
sc3 <- list(s1 = beta_advantage, s2 = 0., s3 = 0.4)
sc4 <- list(s1 = beta_advantage, s2 = 0., s3 = -0.4)
sc5 <- list(s1 = beta_advantage, s2 = 0., s3 = 0)
cols <- c(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:4], mygray)

pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/simulation_growth_WT_M_witherror_all", ".pdf"), width = 4*3*0.8, height = 4*0.8)

par(mfrow = c(1,3), mar = c(4,4,1,1))

for(Nsamplesize in c(30000, 10000, 3000, 30000)){

  
  plot(NULL, xlim = c(-0.25, 0.25), ylim = c(-0.25, 0.25), xlab = "Growth rate historical strains", ylab = "Growth rate emerging variant", las = 1, main = paste0("Daily sample size, N = ", Nsamplesize))
  
  # test a series of (R, duration) combinations that give the same difference in r betwen WT and M (eqaul to 0.07)
  mean_rW <- new_Rtor(param = c(MUW^2/SDW^2, SDW^2/MUW), R = mean(R0vecH)) # mean r for WT
  #target_rM <- mean_rW + 0.07 # target for the delta of growth rates between mutant and WT
  ta_seq <- c()
  
  for(k in 4:1){
    
      myscenario <- get(paste0("sc", k))
      transmission_advantage <- myscenario$s1
      duration_advantage <- myscenario$s2
      sd_advantage <- myscenario$s3
        
      # previously: new_rtoR(param = c(MUM^2/SDM^2, SDM^2/MUM), r = target_rM) / mean(R0vecH) - 1 # transmission advantage needed to reach target rM
      ta_seq <- c(ta_seq, transmission_advantage)
      #plot(siW[1:20], type = "o", pch = 20, ylab = "Distribution serial interval")
      #points(siM[1:20], type = "o", pch = 20, col = "red") 
      
      # need to re-source parameters
      # define distribution of generation time for the mutant by reloading the sourcing file (which updates variant distribution with new values of duration advantage)
      source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")
      
      # simulate:
      sim <- simulate_Rcpp(
        max_age = 60, # maximum age of infection to consider
        initial_IW = 4000, # initial number of infected spread on each age (around 8000 infected per day in IDF if half are detected)
        initial_IM = 400,
        tmax = max_time, # maximum time (3 months)
        simtlockdown = 50, # date of the lockdown (only when step function)
        ifr = 0.0, # ifr
        total_pop = pop_idf * (1 - immunity_idf), # for population immmunity - size of the population left to infect
        par = c(
          R0vecH,
          0.5, # pdetection min
          0.5, # pdetection max
          1,   # k time coefficient
          100, # time mid-point for detection
          transmission_advantage # transmission advantage
        ),
        betafun = 3, # mode for beta function (smooth sigmoid or step function)
        siW = siW, # generation time ("serial interval" distribution)
        siM = siM, # generation time ("serial interval" distribution)
        tid = tid, # time from infection to death
        tit = tit # time from infection to test
      )
      
      
      # cases per day:
      # plot(sim$all_detected_incidenceW, col = "blue", type = "o", pch = 20, lwd = 4, las = 1, ylim = c(0, 10000), ylab = "Detected cases", xlab = "Time (days)", main = paste("transmission advantage = ", transmission_advantage, ", relative duration = ", duration_advantage))
      # points(sim$all_detected_incidenceM, col = "red", type = "o", pch = 20, lwd = 4)
      # points(sim$all_detected_incidenceW+sim$all_detected_incidenceM, col = "gray", type = "o", pch = 20, lwd = 4)
      
      # create a simulated dataset: (note we remove the first 10 days)
      sim_sub_b <- generate_sim_data(mysim = sim, freq_sample_size = Nsamplesize, weekly = F)
      
      points(sim_sub_b$rH_error, sim_sub_b$rV_error, pch = 20, col = cols[k], cex = 0.7)
      abline(0,1)
      
      ll <- loess(sim_sub_b$rV_error ~ sim_sub_b$rH_error); yfit <- predict(ll, newdata = seq(-0.2, 0.1, 0.01))
      # ll <- lm(rV_error ~ rH_error, data = sim_sub_b); yfit <- predict(ll, newdata = data.frame(rH_error = seq(-0.2, 0.1, 0.01)))
      
      lines(seq(-0.2, 0.1, 0.01), yfit, col=cols[k],lty = 1, lwd=2)
      #abline(lm0$coefficients[1], lm0$coefficients[2], col = cols[col_idx])

      
  }
  
  legend("bottomright", col = cols, pch = 20, legend = c("Shorter mgt", "Longer mgt", "Larger var", "Shorter var"), bty = "n")
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)
  if(Nsamplesize==3000) dev.off()
}




# check and visualise the correlation between true R and predicted R DONE

# note the difference in absolute values can be explained mostly by the  the impact of detected individuals and secondarily by the build-up of population immunity
pdf("~/ownCloud/coronavirus/variant_France/clean_code/R0_dynamics_correlation_example.pdf", width = 6, height = 4*2)

# re-simulate with random parameter variation
R0vecH <- c(2.1307553, 2.0105587, 2.0023232, 1.9215775, 1.9642554, 1.9238395, 2.1029669, 2.0780443, 2.1581173, 2.0758108, 2.1253340, 1.9948015, 1.8292173, 1.7309693, 1.7447241, 1.7712299, 1.7653672, 1.7129285, 1.8179133, 1.7841625, 1.9068882, 1.9366379, 2.0484530, 2.1590289, 2.1527725, 2.2736368, 2.2387613, 2.2257805, 2.1303752, 2.1895513, 2.2356056, 2.2355079, 2.3258214, 2.2258386, 2.1631581, 2.2186670, 2.1714491, 2.2069247, 2.1436307, 2.1366790, 2.0831706, 1.9912173, 1.8703819, 1.6309556, 1.4634919, 1.3115412, 1.1725386, 1.0649104, 1.0250046, 0.9133384, 1.0018016, 0.9228869, 0.9090948, 0.7077709, 0.5578058, 0.6371114, 0.7721366, 0.7062939, 0.6210061, 0.4607542, 0.5308621, 0.3945342, 0.4594047, 0.6524794, 0.6660125, 0.5297658, 0.5369410, 0.5645049, 0.6640699, 0.6662376, 0.5892468, 0.5752508, 0.7023259, 0.6168691, 0.4811847, 0.3900102, 0.4270968, 0.2734790, 0.3581384, 0.3733780)
transmission_advantage <- 0.1
duration_advantage <- 0
sd_advantage <- 0
# need to re-source parameters
# define distribution of generation time for the mutant by reloading the sourcing file (which updates variant distribution with new values of duration advantage)
source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")

# simulate:
sim <- simulate_Rcpp(
  max_age = 60, # maximum age of infection to consider
  initial_IW = 4000, # initial number of infected spread on each age (around 8000 infected per day in IDF if half are detected)
  initial_IM = 400,
  tmax = max_time, # maximum time (3 months)
  simtlockdown = 50, # date of the lockdown (only when step function)
  ifr = 0.0, # ifr
  total_pop = pop_idf * (1 - immunity_idf), # for population immmunity - size of the population left to infect
  par = c(
    R0vecH,
    0.5, # pdetection min
    0.5, # pdetection max
    1,   # k time coefficient
    100, # time mid-point for detection
    transmission_advantage # transmission advantage
  ),
  betafun = 3, # mode for beta function (smooth sigmoid or step function)
  siW = siW, # generation time ("serial interval" distribution)
  siM = siM, # generation time ("serial interval" distribution)
  tid = tid, # time from infection to death
  tit = tit # time from infection to test
)
sim_sub_b <- generate_sim_data(mysim = sim, freq_sample_size = Nsamplesize, weekly = F)
# re-simulate with random R0 variation
par(mfrow = c(2,1), mar = c(4,4,1,1))

# cases per day:
max_incidence <- max(sim$all_detected_incidenceW+sim$all_detected_incidenceM)
plot(sim$all_detected_incidenceW, col = "black", type = "l", pch = 20, lwd = 1, lty = 1, las = 1, ylim = c(0, 1.1*max_incidence), ylab = "Detected cases", xlab = "Time (days)", las = 1)
     #,main = paste("transmission advantage = ", signif(transmission_advantage, 3), ", relative duration = ", signif(duration_advantage, 3)))
points(sim$all_detected_incidenceM, col = cols[k], type = "l", lwd = 1, lty = 1)
#points(sim$all_detected_incidenceW+sim$all_detected_incidenceM, col = "gray", type = "l", lwd = 1)
legend("topleft", col = c("black", cols[k]), legend = c("Historical strains cases", "Emerging variant cases"), lty = 2, lwd = 1, bty = "n")

plot(rtoR(r = sim_sub_b$rH, mu = MUW, SD = SDW), type = "l", ylim = c(0., 3.2), xlab = "Time (days)", ylab = expression(paste("Reproduction number ", R)) , lwd = 1, col = "black", las = 1)
points(rtoR(r = sim_sub_b$rV, mu = MUM, SD = SDM), type = "l", lwd = 1, col = cols[k])
points(R0vecH[11:max_time], type = "l", col = "black", lwd = 3)
points(R0vecH[11:max_time] * (1 + transmission_advantage), type = "l", col = cols[k], lwd = 3, lty = 1)
abline(h = 1, lty = 2)
legend("topright", col = c("black", cols[k], "black", cols[k]), lwd = c(3,3,1,1), lty =1 , legend = c(expression(paste("True ", R[0], " historical strains")), expression(paste("True ", R[0], " emerging variant")), expression(paste(R, " realised in H cases")), expression(paste(R, " realised in E cases"))), bty = "n")

# covariance between the two
if(F){
  cc <- ccf(x = R0vecH[11:max_time], y = rtoR(r = sim_sub_b$rH_error, mu = MUW, SD = SDW), plot = F)
  plot(cc$lag, cc$acf, type = "l", ylim = c(0,1), ylab = "Autocorrelation between true and realised Re", xlab = "Lag")
  abline(v = cc$lag[which.max(cc$acf)], lty = 2)
}

dev.off()


################################################## CHECK ANALYTICAL FORMULA FOR COVARIANCE ##################################################


samp_error <- list()
nsamp_total <- 10000
for(i in 1:nsamp_total) samp_error[[i]] <- data.frame(generate_sim_data(mysim = sim, freq_sample_size = Nsamplesize, weekly = F))
idx <- 15 # date index

# useful quantities
I1 <- samp_error[[1]][idx, "Psmoothed"]
I0 <- samp_error[[1]][idx - 1, "Psmoothed"]
p1 <- samp_error[[1]][idx, "fVOC"]
p0 <- samp_error[[1]][idx - 1, "fVOC"]
N1 <- Nsamplesize
N0 <- Nsamplesize

# convergence towards true mean

xseq <- seq(1, nsamp_total, 100)
plot(xseq, sapply(xseq, function(nsamp) mean(unlist(lapply(samp_error[1:nsamp], function(tab) tab[idx, "rH_error"])))), type = "l", ylab = "rH")
abline(h = samp_error[[1]][idx, "rH"], lwd =3, col = "darkgreen")

plot(xseq, sapply(xseq, function(nsamp) mean(unlist(lapply(samp_error[1:nsamp], function(tab) tab[idx, "rV_error"])))), type = "l", ylab = "rV")
abline(h = samp_error[[1]][idx, "rV"], lwd =3, col = "darkgreen")

# variance of the rH, covariance rH rV same time
mean_rH0 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rH_error"])))
mean_rV0 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rV_error"])))
mean_rH1 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx, "rH_error"])))
mean_rV1 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx, "rV_error"])))
mean_rH_squared <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx, "rH_error"]^2)))
mean_rV_squared <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx, "rV_error"]^2)))
mean_rH_rV <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx, "rH_error"] * tab[idx, "rV_error"])))
mean_rH0_rH1 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rH_error"] * tab[idx, "rH_error"])))
mean_rV0_rV1 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rV_error"] * tab[idx, "rV_error"])))
mean_rH0_rV1 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rH_error"] * tab[idx, "rV_error"])))
mean_rH1_rV0 <- mean(unlist(lapply(samp_error[1:nsamp_total], function(tab) tab[idx-1, "rV_error"] * tab[idx, "rH_error"])))

var_rH <- mean_rH_squared - mean_rH1^2
var_rV <- mean_rV_squared - mean_rV1^2
cov_rH_rV <- mean_rH_rV - mean_rH1*mean_rV1
cov_rH0_rH1 <- mean_rH0_rH1 - mean_rH0*mean_rH1
cov_rV0_rV1 <- mean_rV0_rV1 - mean_rV0*mean_rV1
cov_rH0_rV1 <- mean_rH0_rV1 - mean_rH0 * mean_rV1
cov_rH1_rV0 <- mean_rH1_rV0 - mean_rH1 * mean_rV0

pred <- c(
  1/I1 + 1/I0 + p0 / (N0 * (1-p0)) + p1 / (N1 * (1-p1)),  # theoretical prediction var(rH)
  1/I1 + 1/I0 + (1-p0) / (N0 * p0) + (1-p1) / (N1 * p1),  # theoretical prediction var(rV)
  1/I1 + 1/I0 - 1/N0 - 1/N1,                              # theoretical prediction cov(rH, rV)
  -1/I0 - p0/(N0 * (1-p0)),                               # theoretical prediction cov(rH0, rH1)
  -1/I0 - (1-p0)/(N0 * p0),                               # theoretical prediction cov(rV0, rV1)
  -1/I0 + 1/N0,                                           # theoretical prediction cov(rH0, rV1) and  cov(rV0, rH1)
  -1/I0 + 1/N0
)

pdf("~/ownCloud/coronavirus/variant_France/clean_code/check_covariances.pdf", width = 4, height = 4)
par(mar = c(4,4,1,1))
plot(c(var_rH, var_rV, cov_rH_rV, cov_rH0_rH1, cov_rV0_rV1, cov_rH0_rV1, cov_rH1_rV0), pred, pch = 20, xlab = "Simulated variances and covariances", ylab = "Predicted variances and covariances", las = 1)
abline(0, 1, lty = 2)
dev.off()

######################################################################################################################################################
##################################################            INFER FROM SIMULATED DATA             ##################################################
######################################################################################################################################################

initial_da_vec <- seq(-0.4, 0.4, 0.2) 
n_iter_optim <- 20

burnin <- 10
init_IM <- 80
baseline_transmissibility <- 0.3; baseline_transmissibility2 <- 0.45
nrepR0 <- 5

get_initRW <- function(sim_tab, varcov_mat, initial_ta = NULL, initial_da = 0, weekly){
  
  # A function to generate a good starting point for RW
  # Create a good initial RW vector
  # Based on rH and rV and a priori idea of transmission advantage and duration advantage
  if(!weekly){
    weekly_factor <- 1 
    nrows <- max_time
  } else {
    weekly_factor <- 7
    nrows <- nweeks + burnin # the burnin period was already included in nweeks
  }
  initRW1 <- rtoR(sim_tab$rH_error[2:(nrows-burnin)], MUW/weekly_factor, SDW/weekly_factor) # initial values for RW (using data)
  initRW2 <-  rtoR(sim_tab$rV_error[2:(nrows-burnin)], MUW * (1 + initial_da)/weekly_factor, SDW/weekly_factor)
  
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

# npar1 <- 3 + ndays - 1
par(mfrow = c(3,1))
for(sample_size in c(1000, 10000)){

  res <- matrix(NA, nrow = 10*nrepR0*10, ncol = 6 + max_time)
  idx <- 1

  for(k in 1:10){
    
    if(k == 1){ transmission_advantage <- baseline_transmissibility ; duration_advantage <- -0.4 }
    if(k == 2){ transmission_advantage <- baseline_transmissibility ; duration_advantage <- -0.2 }
    if(k == 3){ transmission_advantage <- baseline_transmissibility ; duration_advantage <-  0 }
    if(k == 4){ transmission_advantage <- baseline_transmissibility2 ; duration_advantage <-  0.2 }
    if(k == 5){ transmission_advantage <- baseline_transmissibility2 ; duration_advantage <-  0.4 }
    if(k == 6){ transmission_advantage <- 0.1 ; duration_advantage <- 0 }
    if(k == 7){ transmission_advantage <- 0.2 ; duration_advantage <- 0 }
    if(k == 8){ transmission_advantage <- 0.3 ; duration_advantage <- 0 }
    if(k == 9){ transmission_advantage <- 0.4 ; duration_advantage <- 0 }
    if(k == 10){ transmission_advantage <- 0.5 ; duration_advantage <-0 }
    
    for(repR0 in 1:nrepR0){ # R0 trajectories with different variability (linear decline)
      
      # draw random R0H trajectory:
      #flag <- T
      #cond_on_R0 <- function(myR0vec) mean(myR0vec[1:20]) > 1.1 & mean(myR0vec[(max_time-20):max_time]) < 0.9 & all(myR0vec > 0.5) & all(myR0vec < 5)
      #while(flag){
      #  R0vecH <- mvrnorm(n = 1, mu = rep(R0H, max_time), Sigma = create_BMautocov_mat(temporal_autocor, max_time = max_time))
      #  if(cond_on_R0(R0vecH)) flag <- F # select a realisation of BM where it is below 1 at the end
      #}
      initial_R0 <- 1 + repR0 / 10
      final_R0 <- 1 - repR0 / 10
      R0vecH <- c(rep(initial_R0, 10), seq(initial_R0, final_R0, by = -(initial_R0 - final_R0)/(70-1)))
      
      cat(transmission_advantage, duration_advantage, "\n")
      #print(R0vecH)
      # define distribution of generation time for the mutant by reloading the sourcing file (which updates variant distribution with new values of duration advantage)
      source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")
      
      # simulate:
      sim <- simulate_Rcpp(
        max_age = 60, # maximum age of infection to consider
        initial_IW = 4000, # initial number of infected spread on each age (around 8000 infected per day in IDF if half are detected)
        initial_IM = init_IM,
        tmax = max_time, # maximum time (3 months)
        simtlockdown = 50, # date of the lockdown (only when step function)
        ifr = 0.0, # ifr
        total_pop = pop_idf * (1 - immunity_idf), # for population immmunity - size of the population left to infect
        par = c(
          R0vecH,
          0.5, # pdetection min
          0.5, # pdetection max
          1,   # k time coefficient
          100, # time mid-point for detection
          transmission_advantage # transmission advantage
        ),
        betafun = 3, # mode for beta function (smooth sigmoid or step function)
        siW = siW, # generation time ("serial interval" distribution)
        siM = siM, # generation time ("serial interval" distribution)
        tid = tid, # time from infection to death
        tit = tit # time from infection to test
      )
      # note: incidence infections at time i (all_incidenceW and all_incidenceM) reflect transmission that occurred between i-1 and i
      # the first number of all_detected_incidenceM is initial_IM *  pdetection because of the assumption that initially all ages of infection are equally representd
      # cases are shifted compared to infections
      
      for(m in 1:10){ # replicates
        
        # create the simulated dataset from the simulations (with error)
        sim_sub_b <- generate_sim_data(mysim = sim, freq_sample_size = sample_size, weekly = F)
        names(sim_sub_b)[names(sim_sub_b) == "N"] <- "Ntrue"
        if(any(is.infinite(sim_sub_b$rH_error)) | any(is.infinite(sim_sub_b$rV_error))) next # go to next one if infinite rH or rV
        max_frequency_reached <- max(sim_sub_b$fVOC, na.rm = T)
        
        # Var-covar matrix of the multivariate normal distribution
        varcov <- create_cov_mat(tab = sim_sub_b)
        stopifnot(isSymmetric.matrix(varcov)) #  check matrix is symmetric:
        
        # data
        mydata <- c(rbind(sim_sub_b$rH_error[2:(max_time-burnin)], sim_sub_b$rV_error[2:(max_time-burnin)])) # estimated r; what we are trying to fit
        
        # the true_RW is:
        true_RW1 <- rtoR(sim_sub_b$rH[2:(max_time-burnin)], MUW, SDW) 
        true_RW2 <- rtoR(sim_sub_b$rV[2:(max_time-burnin)], MUM, SDM) / (1 + transmission_advantage)
        true_RW <- 0.5 * (true_RW1 + true_RW2)
        
        
        opt2_list <- list()
        RW_error <- c()
        idx_list <- 1
        
        for(my_initial_da in initial_da_vec){
          
          # generate a good initRW:
          tmp <- get_initRW(sim_tab = sim_sub_b, varcov_mat = varcov, initial_ta = NULL, initial_da = my_initial_da, weekly = FALSE)
          initRW <- (tmp$initRW)
          initial_ta <- tmp$initial_ta
          initial_da <- tmp$initial_da
          if(any(is.na(initRW))) stop("Some elements of initRW are NA")
        
          lower2 <- c(0, -0.5)
          upper2 <- c(1, 1)
          npar2 <- 2
          
          initfun2 <- function() c(initial_ta, initial_da)
          opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = myneglik2, init.fun = initfun2, lower = lower2, upper = upper2, mydata = mydata, myinitRW = initRW, vars = varcov, fixed_s3 = 0, weekly = F)
          opt2_list[[idx_list]] <- opt2; idx_list <- idx_list + 1
          RW_error <- c(RW_error, mean((initRW - true_RW)^2))
          #2nd round:
          
          for(l in 1:n_iter_optim){
            tmp <- get_initRW(sim_tab = sim_sub_b, varcov_mat = varcov, initial_ta = opt2$pars[1], initial_da = opt2$pars[2], weekly = FALSE)
            initRW <- (tmp$initRW)
            initial_ta <- tmp$initial_ta
            initial_da <- tmp$initial_da
            if(initial_da >= upper2[2]) initial_da <- upper2[2] * 0.9
            if(initial_ta >= upper2[1]) initial_ta <- upper2[1] * 0.9
            initfun2 <- function() c(initial_ta, initial_da)
            opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = myneglik2, init.fun = initfun2, lower = lower2, upper = upper2, mydata = mydata, myinitRW = initRW, vars = varcov, fixed_s3 = 0, weekly = F)
            opt2_list[[idx_list]] <- opt2; idx_list <- idx_list+1
            RW_error <- c(RW_error, mean((initRW - true_RW)^2))
          }
        }
        
        all_lik <- unlist(lapply(opt2_list, function(x) x$lik))
        opt2 <- opt2_list[[which.min(all_lik)]]
        
        # plot a few things
        plot(all_lik, type = "o", pch = 20,  main = paste("max freq", signif(max_frequency_reached, 3), " - best par:", paste(signif(opt2$pars, 3), collapse = ",  ")))
        plot(initRW, type = "l", lwd = 3, pch= 20, col= "blue", ylim = c(0,2))
        points(true_RW, type = "l", col = "gray", lwd = 3)
        plot(all_lik, RW_error, pch = 20)
        points(all_lik[which.min(all_lik)], RW_error[which.min(all_lik)], pch = 20, cex = 4, col = "red")
        # save results
        res[idx, 1:6] <- c(transmission_advantage, duration_advantage, opt2$pars, opt2$lik, max_frequency_reached)
        res[idx, 7:(7+max_time-1)] <- R0vecH
        idx <- idx + 1
        
      } # end of loop on replicates with error
    } # end of loop on R0 trajectory more or less variable
  } # end of loop on parameters

  res <- data.frame(res)
  names(res) <- c("ta", "da", "inf_ta", "inf_da", "max_lik", "max_fVOC", paste0("R0", 1:max_time))
  # plot(res$ta, res$inf_ta); abline(0,1)
  # plot(res$da, res$inf_da); abline(0,1)
  
  assign(x = paste0("res_", sample_size), res)

}



######################################################################################################################################################
######################################            INFER FROM WEEKLY SIMULATED DATA             #######################################################
######################################################################################################################################################

burnin <- 1
sample_size <- 1000

res <- matrix(NA, nrow = 10*nrepR0*10, ncol = 6 + max_time)
idx <- 1

for(k in 1:10){
  
  if(k == 1){ transmission_advantage <- baseline_transmissibility ; duration_advantage <- -0.4 }
  if(k == 2){ transmission_advantage <- baseline_transmissibility ; duration_advantage <- -0.2 }
  if(k == 3){ transmission_advantage <- baseline_transmissibility ; duration_advantage <-  0 }
  if(k == 4){ transmission_advantage <- baseline_transmissibility2 ; duration_advantage <-  0.2 }
  if(k == 5){ transmission_advantage <- baseline_transmissibility2 ; duration_advantage <-  0.4 }
  if(k == 6){ transmission_advantage <- 0.1 ; duration_advantage <- 0 }
  if(k == 7){ transmission_advantage <- 0.2 ; duration_advantage <- 0 }
  if(k == 8){ transmission_advantage <- 0.3 ; duration_advantage <- 0 }
  if(k == 9){ transmission_advantage <- 0.4 ; duration_advantage <- 0 }
  if(k == 10){ transmission_advantage <- 0.5 ; duration_advantage <-0 }
  
  for(repR0 in 1:nrepR0){ # R0 trajectories with different variability (linear decline)
    
    # draw random R0H trajectory:
    #flag <- T
    #while(flag){
    #  R0vecH <- mvrnorm(n = 1, mu = rep(R0H, max_time), Sigma = create_BMautocov_mat(temporal_autocor, max_time = max_time))
    #  if(cond_on_R0(R0vecH)) flag <- F # select a realisation of BM where it is below 1 at the end
    #}
    initial_R0 <- 1 + repR0 / 10
    final_R0 <- 1 - repR0 / 10
    R0vecH <- c(rep(initial_R0, 10), seq(initial_R0, final_R0, by = -(initial_R0 - final_R0)/(70-1)))
    cat(transmission_advantage, duration_advantage, "\n")
    #print(R0vecH)
    
    # define distribution of generation time for the mutant by reloading the sourcing file (which updates variant distribution with new values of duration advantage)
    source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")
    
    # simulate:
    sim <- simulate_Rcpp(
      max_age = 60, # maximum age of infection to consider
      initial_IW = 4000, # initial number of infected spread on each age (around 8000 infected per day in IDF if half are detected)
      initial_IM = init_IM,
      tmax = max_time, # maximum time (3 months)
      simtlockdown = 50, # date of the lockdown (only when step function)
      ifr = 0.0, # ifr
      total_pop = pop_idf * (1 - immunity_idf), # for population immmunity - size of the population left to infect
      par = c(
        R0vecH,
        0.5, # pdetection min
        0.5, # pdetection max
        1,   # k time coefficient
        100, # time mid-point for detection
        transmission_advantage # transmission advantage
      ),
      betafun = 3, # mode for beta function (smooth sigmoid or step function)
      siW = siW, # generation time ("serial interval" distribution)
      siM = siM, # generation time ("serial interval" distribution)
      tid = tid, # time from infection to death
      tit = tit # time from infection to test
    )
    
    for(m in 1:10){ # replicates on data
      
      # create the simulated dataset from the simulations (with error)
      sim_sub_b <- generate_sim_data(mysim = sim, freq_sample_size = sample_size, weekly = T)
      names(sim_sub_b)[names(sim_sub_b) == "N"] <- "Ntrue"
      if(any(is.infinite(sim_sub_b$rH_error)) | any(is.infinite(sim_sub_b$rV_error))) next # go to next one if infinite rH or rV
      nweeks <- nrow(sim_sub_b)
      max_frequency_reached <- max(sim_sub_b$fVOC, na.rm = T)
      
      # Var-covar matrix of the multivariate normal distribution
      varcov_weekly <- create_cov_mat(tab = sim_sub_b)
      stopifnot(isSymmetric.matrix(varcov_weekly))
      
      # data
      mydata <- c(rbind(sim_sub_b$rH_error[2:nweeks], sim_sub_b$rV_error[2:nweeks])) # estimated r; what we are trying to fit
      
      opt2_list <- list()
      RW_error <- c()
      idx_list <- 1
      
      for(my_initial_da in initial_da_vec){
        
        # generate a good initRW:
        tmp <- get_initRW(sim_tab = sim_sub_b, varcov_mat = varcov_weekly, initial_da = my_initial_da, weekly = TRUE)
        initRW <- tmp$initRW
        initial_ta <- tmp$initial_ta
        initial_da <- tmp$initial_da
        
        if(any(is.na(initRW))) stop("Some elements of initRW are NA")
        
        # the true_RW is:
        true_RW1 <- rtoR(sim_sub_b$rH[2:nweeks], MUW/7, SDW/7) 
        true_RW2 <- rtoR(sim_sub_b$rV[2:nweeks], MUM/7, SDM/7) / (1 + transmission_advantage)
        true_RW2[is.na(true_RW2)] <- true_RW1[is.na(true_RW2)]
        true_RW <- 0.5 * (true_RW1 + true_RW2)
        
        if(any(is.na(initRW))) stop("Some elements of initRW are NA")
        
        #plot(initRW, type = "l", lwd = 3, pch= 20, col= "blue", ylim = c(0,5))
        #points(true_RW, type = "l", col = "gray", lwd = 3)
      
        initfun2 <- function() c(initial_ta, initial_da)
        opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = myneglik2, init.fun = initfun2, lower = lower2, upper = upper2, mydata = mydata, myinitRW = initRW, vars = varcov_weekly, fixed_s3= 0, weekly = T)
        
        # 2nd round
        opt2_list[[idx_list]] <- opt2; idx_list <- idx_list + 1
        RW_error <- c(RW_error, mean((initRW - true_RW)^2))
        
        for(l in 1:n_iter_optim){
          tmp <- get_initRW(sim_tab = sim_sub_b, varcov_mat = varcov_weekly, initial_ta = opt2$pars[1], initial_da = opt2$pars[2], weekly = TRUE)
          initRW <- tmp$initRW
          initial_ta <- tmp$initial_ta
          initial_da <- tmp$initial_da
          initfun2 <- function() c(initial_ta, initial_da)
          opt2 <- optim.fun.repeated(n.repeats = 1, lik.fun = myneglik2, init.fun = initfun2, lower = lower2, upper = upper2, mydata = mydata, myinitRW = initRW, vars = varcov_weekly,  fixed_s3= 0, weekly = T)
          opt2_list[[idx_list]] <- opt2; idx_list <- idx_list + 1
          RW_error <- c(RW_error, mean((initRW - true_RW)^2))
        }
      }
      all_lik <- unlist(lapply(opt2_list, function(x) x$lik))
      opt2 <- opt2_list[[which.min(all_lik)]]
      res[idx, 1:6] <- c(transmission_advantage, duration_advantage, opt2$pars, opt2$lik, max_frequency_reached)
      res[idx, 7:(7+max_time-1)] <- R0vecH
      idx <- idx + 1
      
      # ad plot a few things
      plot(all_lik, type = "o", pch = 20,  main = paste("max freq", signif(max_frequency_reached, 3), " - best par:", paste(signif(opt2$pars, 3), collapse = ",  ")))
      plot(initRW, type = "l", lwd = 3, pch= 20, col= "blue", ylim = c(0,2))
      points(true_RW, type = "l", col = "gray", lwd = 3)
      plot(all_lik, RW_error, pch = 20)
      points(all_lik[which.min(all_lik)], RW_error[which.min(all_lik)], pch = 20, cex = 4, col = "red")
      
    } # end of loop on replicates with error
  } # end of loop on R0 trajectories
} # end of loop on parameters s1, s2

res <- data.frame(res)
names(res) <- c("ta", "da", "inf_ta", "inf_da", "max_lik", "max_fVOC", paste0("R0", 1:max_time))
assign(x = paste0("res_", sample_size, "_weekly"), res)

save.image(file = "~/ownCloud/coronavirus/variant_France/clean_code/simulations_30August2021.RData")




##################################################                    PLOT THE OVERALL FIGURE                  ##################################################

# load("~/ownCloud/coronavirus/variant_France/clean_code/simulations_30August2021.RData")

cols2 <- RColorBrewer::brewer.pal(n = 3, name = "Set1")

pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/simulation_study_inference.pdf"), width = 4*3*0.6, height = 3*4.2*0.6)

par(mfrow =c(3,2), mar = c(5,4,2,1))
for(initial_R0 in c(1.5, 1.3, 1.1)){
  ylims <- c(0, 0.6)
  plot(NULL, xlim = ylims, ylim  = ylims,
       xlab = "True transmission advantage", ylab = "Inferred transmission advantage", las = 1, bty = "n")
  abline(0, 1, col = "gray")
  idx <- 1
  
  for(tab in c("res_1000", "res_10000", "res_1000_weekly")){
    truc <- get(tab)
    ylims <- c(0, 0.6)
    sub <- which(truc$R01==initial_R0 & truc$da == 0)
    points(truc$ta[sub] + (idx-1)/100., truc$inf_ta[sub], pch = 20, col = cols2[idx])
    idx <- idx + 1
  }
  legend("topleft", legend = c("1000 daily", "10000 daily", "7000 weekly"), col = cols2, pch = 20, bty = "n")
  
  ylims <- c(-0.6, 1)
  plot(NULL, xlim = ylims, ylim  = ylims,
       xlab = "True relative mgt", ylab = "Inferred relative mgt", las = 1, bty = "n")
  abline(0, 1, col = "gray")
  idx <- 1
  for(tab in c("res_1000", "res_10000", "res_1000_weekly")){
    truc <- get(tab)
    sub <- which(truc$R01==initial_R0 & (truc$ta == baseline_transmissibility | truc$ta == baseline_transmissibility2))
    points(truc$da[sub] + 2*(idx-1)/100., truc$inf_da[sub], pch = 20, col = cols2[idx])
    idx <- idx + 1
  }
  legend("topleft", legend = c("1000 daily", "10000 daily", "7000 weekly"), col = cols2, pch = 20, bty = "n")
}

dev.off()





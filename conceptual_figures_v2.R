
figure_size_factor <- 1.4
figure_size_factor2 <- 1.6
mygray <- "gray60"

# must enter the first 18 lines of "inference_UK_France_vX.R" script
# conceptual figures

rtoR(r = 0.05, mu = MUW, SD = SDW)
Rtor(R = 2.5, mu = MUW, SD = SDW)

mylegend <- c("Historical", "Shorter mean generation time", "Longer mean generation time", "Larger variance", "Smaller variance")

for(metascenario in 1:2){
  
  if(metascenario==1){
    sc1 <- list(s1 = 0., s2 = 5.5/6.5-1, s3 = 0.)
    sc2 <- list(s1 = 0., s2 = 7.5/6.5-1, s3 = 0.)
    sc3 <- list(s1 = 0., s2 = 0., s3 = 0.4)
    sc4 <- list(s1 = 0., s2 = 0., s3 = -0.4)
    sc5 <- list(s1 = 0., s2 = 0., s3 = 0)
  }
  if(metascenario==2){
    beta_advantage <- 0.1
    sc1 <- list(s1 = beta_advantage, s2 = 5.5/6.5-1, s3 = 0.)
    sc2 <- list(s1 = beta_advantage, s2 = 7.5/6.5-1, s3 = 0.)
    sc3 <- list(s1 = beta_advantage, s2 = 0., s3 = 0.4)
    sc4 <- list(s1 = beta_advantage, s2 = 0., s3 = -0.4)
    sc5 <- list(s1 = beta_advantage, s2 = 0., s3 = 0)
  }
  cols <- c(RColorBrewer::brewer.pal(n = 8, name = "Paired")[1:4], mygray)
  
  RW <- 1
  
  get_mutant_par <- function(mysc){
    return(
      list(
        RV = with(mysc, RW*(1+s1)),
        MUV = with(mysc, MUW * (1+s2)),
        SDV = with(mysc, SDW * (1+s3)),
        shapeV = with(mysc, (MUW * (1+s2))^2/(SDW * (1+s3))^2),
        scaleV = with(mysc, (SDW * (1+s3))^2/( MUW * (1+s2)))
      )
    )
  }
  
  x <- seq(0, 20, 0.1)
  
  pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figure_gamma_distribution_", metascenario, ".pdf"), width = 4*figure_size_factor, height = 3*figure_size_factor)
  par(mar = c(4,4,1,1))
  plot(x, RW * dgamma(x, shape = MUW^2/SDW^2, scale = SDW^2/MUW), type = "l", xlab ="Days post infection", ylab  = "Density", xaxs = "i", yaxs = "i", las = 1, lwd = 3, ylim = c(0, 0.25), bty = "n")
  with(get_mutant_par(sc1), points(x, RV * dgamma(x, shape = shapeV, scale = scaleV), type = "l", col = cols[1], lwd = 3))
  with(get_mutant_par(sc2), points(x, RV * dgamma(x, shape = shapeV, scale = scaleV), type = "l", col = cols[2], lwd = 3))
  with(get_mutant_par(sc3), points(x, RV * dgamma(x, shape = shapeV, scale = scaleV), type = "l", col = cols[3], lwd = 3))
  with(get_mutant_par(sc4), points(x, RV * dgamma(x, shape = shapeV, scale = scaleV), type = "l", col = cols[4], lwd = 3))
  #if(metascenario==2)with(get_mutant_par(sc5), points(x, RV * dgamma(x, shape = shapeV, scale = scaleV), type = "l", col = cols[5], lwd = 3))
  if(metascenario == 1) legend("topright", legend = mylegend, lwd = 3, col = c("black", cols), bty = "n")
  dev.off()
  
  # figure selection coefficient for the three scenarios:
  
  Rfactor <- seq(0.1, 2, 0.01)
  
  pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figure_r_s_", metascenario, ".pdf"), width = 4*figure_size_factor, height = 1*figure_size_factor*3)
  par(mar = c(4,4,1,1), mfrow=  c(1,1))
  
  #xlab <- expression(paste("Effective reproduction number ", R, " wild type"))
  xlab <- expression(paste("Proxy for transmissibility ", R[H]))
  

  plot(NULL, lty = 1, type = "l", 
         xlab = xlab, ylab = "Selection coefficient",  xaxs = "i", yaxs = "i", las = 1, lwd = 3, ylim= c(-0.05, 0.1), xlim = c(0, 2), axes = F)
  axis(side = 1, at = c(0, 0.5, 1, 1.5, 2))
  axis(side = 2, at = c(-0.04, 0, 0.04, 0.08), las = 1)
  #if(metascenario==2)
  {
    for(mys in c(0, 0.1, 0.2, 0.3, 0.4, 0.5)){
      points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(list(s1 = mys, s2 = 0, s3 = 0)), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = "gray", lty = 1, type = "l", lwd = 1)
    }
  }
  abline(h = 0, lty = 1)
  points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(sc1), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = cols[1], lty = 1, type = "l", lwd = 3)
  points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(sc2), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = cols[2], lty = 1, type = "l", lwd = 3)
  points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(sc3), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = cols[3], lty = 1, type = "l", lwd = 3)
  points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(sc4), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = cols[4], lty = 1, type = "l", lwd = 3)
  #points(Rfactor * RW, y <- sapply(Rfactor, function(f) with(get_mutant_par(sc5), Rtor(R = f*RV, mu = MUV, SD = SDV)) - Rtor(R = f*RW, mu = MUW, SD = SDW)), col = cols[5], lty = 1, type = "l", lwd = 3)
  
  text(x = rep(1.85, 6), y = c(0.003, 0.021, 0.04, 0.057, 0.073, 0.09), labels = c("WT", "+10%", "+20%", "+30%", "+40%", "+50%"), col = mygray)
  
  dev.off()
  
  pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figure_r_r_", metascenario, ".pdf"), width = 4, height = 4)
  
  par(mar = c(4,4,1,1), mfrow = c(1,1))
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
  points(x1[sub1], y1[sub1], col = cols[1], lty = 1, type = "l", lwd = 3)
  points(x2[sub2], y2[sub2], col = cols[2], lty = 1, type = "l", lwd = 3)
  
  points(x = 7 * Rtor(R = 1/(1 + sc1$s1), mu = MUW, SD = SDW), y = 0, pch = 20, cex = 2)
  #text(x = rep(1.85, 6), y = c(0.003, 0.021, 0.04, 0.057, 0.073, 0.09), labels = c("WT", "+10%", "+20%", "+30%", "+40%", "+50%"), col = mygray)
  legend("bottomright", legend = c("Shorter mgt", "Longer mgt"), lwd = 3, col = c(cols[1:2]), bty = "n")
  
  segments(x0 = x2[100], y0 = y2[100], x1 = x2[120], y1 = y2[100], lwd = 1.5, col = cols[2])
  segments(x0 = x2[120], y0 = y2[100], x1 = x2[120], y1 = y2[120], lwd = 1.5, col = cols[2])
  
  text(x = -0.2, y = 0.3, labels = expression(paste("Slope", phantom() %~~% phantom(), "1 + ", s[1], cv[H]^2, " - ", s[2])), col = cols[2])
  
  dev.off()
  
  ######################################################################
  #####                 epidemiological dynamics                   #####
  ######################################################################
  
  R0H <- 2
  max_time <- 80 # duration of simulations (we will eliminate ten days of burnin period to eliminate transient effects)
  max_age <- 60 # following infections up to age max_age
  begin_time <- 1
  
  flag <- F; idx <- 1
  while(flag){
    R0vecH <- mvrnorm(n = 1, mu = rep(R0H, max_time), Sigma = create_BMautocov_mat(0.1, max_time = max_time))
    if(mean(R0vecH[1:20]) > 1.5 & mean(R0vecH[(max_time-20):max_time]) < 0.95  & all(R0vecH > 0.)) flag <- F # select a realisation of BM where it is below 1 at the end
    print(idx <- idx+1)
  }
  R0vecH <- c(2.1307553, 2.0105587, 2.0023232, 1.9215775, 1.9642554, 1.9238395, 2.1029669, 2.0780443, 2.1581173, 2.0758108, 2.1253340, 1.9948015, 1.8292173, 1.7309693, 1.7447241, 1.7712299, 1.7653672, 1.7129285, 1.8179133, 1.7841625, 1.9068882, 1.9366379, 2.0484530, 2.1590289, 2.1527725, 2.2736368, 2.2387613, 2.2257805, 2.1303752, 2.1895513, 2.2356056, 2.2355079, 2.3258214, 2.2258386, 2.1631581, 2.2186670, 2.1714491, 2.2069247, 2.1436307, 2.1366790, 2.0831706, 1.9912173, 1.8703819, 1.6309556, 1.4634919, 1.3115412, 1.1725386, 1.0649104, 1.0250046, 0.9133384, 1.0018016, 0.9228869, 0.9090948, 0.7077709, 0.5578058, 0.6371114, 0.7721366, 0.7062939, 0.6210061, 0.4607542, 0.5308621, 0.3945342, 0.4594047, 0.6524794, 0.6660125, 0.5297658, 0.5369410, 0.5645049, 0.6640699, 0.6662376, 0.5892468, 0.5752508, 0.7023259, 0.6168691, 0.4811847, 0.3900102, 0.4270968, 0.2734790, 0.3581384, 0.3733780)
  plot(R0vecH, type = "l", las= 1, axes = F, xlab = "Time", ylab = expression(paste(R[0]))); abline(h = 1, lty = 2)
  axis(side = 1, at = seq(0, 80, 20))
  axis(side = 2, at = c(0.5, 1, 1.5, 2), las = 1)
  
  # parameterise for IDF
  pop_idf <- 12278210 # pop IDF
  immunity_idf <- 0.2549 # estimation of immunity in IDF on 2021-01-01
  sourceCpp("~/ownCloud/coronavirus/variant_France/clean_code/simulate_v2.cpp")  # simulations in Rcpp
  
  # tune parameters to value
  
  pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figure_evo_", metascenario, ".pdf"), width = 4*figure_size_factor2, height = 3*figure_size_factor2)
  par(mar = c(4,4,1,1), mfrow = c(1,1))
  plot(NULL, xlim = c(0, max_time), ylim = c(0, 1.), ylab = "Frequency of variants", xlab = "Time (days)", las = 1, xaxs = "i", yaxs = "i", axes = F)
  
  # add the simulations with some transmissibility advtange:
  ta_labels <- c("+10%", "+20%", "+30%", "+40%", "+50%")
  names(ta_labels) <- seq(0.1, 0.5, 0.1)
  
  # first the 
  for(transmission_advantage in  seq(0.1, 0.5, 0.1)){
    
    duration_advantage <- 0
    sd_advantage <- 0
    mycol <- "gray"
    source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")
    
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
    f <- sim$all_detected_incidenceM/(sim$all_detected_incidenceW+sim$all_detected_incidenceM)
    points(begin_time:max_time, f[begin_time:max_time], col = mycol, type = "l", lty = 1, pch = 20, lwd = 1)
    text(75, (sim$all_detected_incidenceM/(sim$all_detected_incidenceW+sim$all_detected_incidenceM))[75], labels = ta_labels[as.character(transmission_advantage)], col = mygray)
  }
  all_epi_curves <- list()
  for(scenar in 1:4){
    
    transmission_advantage <- get(paste0("sc", scenar))$s1
    duration_advantage <- get(paste0("sc", scenar))$s2
    sd_advantage <- get(paste0("sc", scenar))$s3
    mycol <- cols[scenar]
    
    source("~/ownCloud/coronavirus/variant_France/clean_code/source_simulation_parameters.R")
    
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
    # plot frequency dynamics:
    f <- sim$all_detected_incidenceM/(sim$all_detected_incidenceW+sim$all_detected_incidenceM)
    points(begin_time:max_time, f[begin_time:max_time], col = mycol, type = "l", pch = 20, lwd = 3, lty = 1) # frequency
    #if(scenar == 1) legend("topleft", col = c("black", cols), legend = c("Wild Type cases", paste0("Variant cases - ", mylegend[2:5])), lty = 1, lwd = 2, bty = "n", cex = 0.5)
    
    all_epi_curves[[scenar]] <- sim
  }
  axis(side = 1, at = seq(0, 80, 40))
  axis(side = 2, at = seq(0, 1, 0.2), las = 1)
  
  # add R0(t)
  old_par <- par()
  par(new = TRUE, mar=c(0,0,2,0), fig = c(0.22,0.48,0.7,1.)) # overlay existing plot
  plot(R0vecH, type = "l", las= 1, lwd = 2, col = "black", col.main = "black", axes = F, xlab = "Time", main = expression(paste(R[0,WT](t))))
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
  points(R0vecH, type = "l", lwd = 2, col = "black")
  axis(side = 2, at = c(0.5, 1, 1.5, 2), las = 1)
  axis(side = 1, at = seq(0, 80, 40), labels = rep("", 3))
  suppressWarnings(par(old_par))
  
  dev.off()
  
  # now to the epidemiological figure
  pdf(paste0("~/ownCloud/coronavirus/variant_France/clean_code/conceptual_figure_epi_", metascenario, ".pdf"), width = 4*figure_size_factor2, height = 3*figure_size_factor2)
  
  par(mar = c(4,5,1,1), mfrow = c(1,1))
  plot(NULL, xlim = c(0, max_time), ylim = c(0, 50000), ylab = "", xlab = "Time (days)", las = 1, xaxs = "i", yaxs = "i", axes = F)
  title(ylab="Number of cases", line=4)
  
  # cases per day:
  for(scenar in 1:4){
    mycol <- cols[scenar]
    points(begin_time:max_time, all_epi_curves[[scenar]]$all_detected_incidenceM[begin_time:max_time] , col = mycol, type = "l", pch = 20, lwd = 3, lty = 1) # mutant incidence
    points(begin_time:max_time, all_epi_curves[[scenar]]$all_detected_incidenceW[begin_time:max_time] , col = "black", type = "l", pch = 20, lwd = 3, lty = 1) # WT incidence
  }
  axis(side = 1, at = seq(0, 80, 40))
  axis(side = 2, at = seq(0, 50000, 10000), las = 1)
  if(metascenario == 1) legend("topleft", legend = mylegend, lwd = 3, col = c("black", cols), bty = "n")
  dev.off()

}




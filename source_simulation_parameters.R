


max_age <- 60 # following infections up to age max_age
max_time <- 80 # duration of simulations (we will eliminate ten days of burnin period to eliminate transient effects)

# 1. discretised serial interval for historical strains
MUW
SDW
siW <- pgamma(1.5:(max_age+0.5), shape = MUW^2/SDW^2, scale = SDW^2/MUW) - pgamma(0.5:(max_age-0.5), shape = MUW^2/SDW^2, scale = SDW^2/MUW); siW[1] <- siW[1] + pgamma(0.5,  shape = MUW^2/SDW^2, scale = SDW^2/MUW) # gamma with mean 6.5, sd 4
stopifnot(abs(sum(siW)-1) < 1e-4)

# 2. discretise the distribution of time from infection to symptom onset
tio <- plnorm(1.5:(max_age+0.5), meanlog = 1.51804, sdlog = 0.471594) - plnorm(0.5:(max_age-0.5), meanlog = 1.51804, sdlog = 0.471594); tio[1] <- tio[1] + plnorm(0.5, meanlog = 1.51804, sdlog = 0.471594) # Lauer et al Annals of Internal Medicine; or Linton et al (meanlog = 1.525, sdlog = 0.6288)
if(sum(tio) > 0.999) tio <- tio / sum(tio) else stop()

# 3. discretise the distribution of time from onset to death
tod <- pgamma(1.5:(max_age+0.5), shape = 5, rate = 1/4) - pgamma(0.5:(max_age-0.5), shape = 5, rate = 1/4); tod[1] <- tod[1] + pgamma(0.5, shape = 5, rate = 1/4)  # Wu et al Nature Medicine
if(sum(tod) > 0.999) tod <- tod / sum(tod) else stop()

# 4 time from infection to death is the convolution of the two
tid <- convolve(tio, rev(tod), type = "o")[1:max_age]
if(sum(tid) > 0.99) tid <- tid / sum(tid) else stop()

# 5 discretise the distribution of time from onset to test and deduce infection to test
tot <- pgamma(1.5:(max_age+0.5), shape = 0.69, rate = 0.31) - pgamma(0.5:(max_age-0.5), shape = 0.69, rate = 0.31); tot[1] <-tot[1] + pgamma(0.5, shape = 0.69, rate = 0.31)  # numerised data from Lauer et al.
sum(tot * 1:max_age)

tit <- convolve(tio, rev(tot), type = "o")[1:max_age]
if(sum(tit) > 0.99) tit <- tit / sum(tit) else stop()



# Finally, define distribution of generation time for the mutant:

MUM <- MUW * (1 + duration_advantage)
SDM <- SDW * (1 + sd_advantage) # fixed to the same value for now
siM <- pgamma(1.5:(max_age+0.5), shape = MUM^2/SDM^2, scale = SDM^2/MUM) - pgamma(0.5:(max_age-0.5), shape = MUM^2/SDM^2, scale = SDM^2/MUM); siM[1] <- siM[1] + pgamma(0.5,  shape = MUM^2/SDM^2, scale = SDM^2/MUM) # gamma with mean 6.5, sd 4
stopifnot(abs(sum(siM)-1) < 1e-4)




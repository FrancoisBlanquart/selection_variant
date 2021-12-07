library(plyr)

##############################                                            DEFINE CONSTANT                                     ############################## 

# TODO CHANGE THIS IF WE CHANGE TIME UNIT (I.E. SHIFT FROM DAY UNIT TO WEEK UNIT)
MUW <- 6.5 # this has unit of time
SDW <- 4 # this has unit of time
temporal_autocor <- 0.05 # this has unit of time

############################################################ UK : MERGE CASES WITH FREQUENCY TABLE AND COMPUTE RH, RV ############################################################

uk_cases <- read.csv("~/ownCloud/coronavirus/variant_France/clean_code/cleaned_UK_cases.csv")
uk_freqs <- read.csv("~/ownCloud/coronavirus/variant_France/clean_code/cleaned_UK_frequencies.csv")
uk_cases$week <- as.character(uk_cases$week)
uk_freqs$week <- as.character(uk_freqs$week)
uk <- merge(uk_cases, uk_freqs, by = c("region", "week"), all = T)
uk <- uk[with(uk, order(region,week)), ]
uk <- uk[!is.na(uk$cases) & !is.na(uk$p), ]
stopifnot(all(uk$day.x==uk$day.y))
stopifnot(all(uk$S_minus+uk$S_plus==uk$N, na.rm =T))
names(uk) <- c("reg2", "jour", "week_number", "weekly_cases", "day", "S_plus", "NVOC", "fVOC", "N", "logitp", "day.y")
uk <- uk[c("reg2", "jour", "week_number", "weekly_cases", "day", "NVOC", "fVOC", "N")]
for(mycol in c("weekly_cases", "NVOC", "N")) uk[mycol] <- as.numeric(unlist(uk[mycol]))
uk_region <- uk
london <- uk_region[uk_region$reg2=="London", ]
uk <- ddply(uk_region, .(week_number), summarise, weekly_cases = sum(weekly_cases), NVOC = sum(NVOC), N = sum(N), jour = unique(jour))

# period 2: growth of the delta variant
uk_freqs2 <- read.csv("~/ownCloud/coronavirus/variant_France/clean_code/cleaned_UK_frequencies_period2.csv")
uk_freqs2$week <- as.character(uk_freqs2$week)
uk_cases2 <- ddply(uk_cases, .(week), summarise, cases = sum(cases), day = unique(day), week_number = unique(week_number))
uk_cases2$week <- as.Date(uk_cases2$week) - 1 # shift one day the cases to match the frequencies (as weekly frequencies start on Monday, cases start on Tuesday)
uk_cases2$day <- uk_cases2$day - 1
uk_cases2$region <- "UK"
uk2 <- merge(uk_cases2, uk_freqs2, by = c("region", "week"), all = T)
uk2 <- uk2[!is.na(uk2$cases) & !is.na(uk2$p), ]
stopifnot(all(uk2$day.x==uk2$day.y))
names(uk2)[names(uk2) == "cases"] <- "weekly_cases"
names(uk2)[names(uk2) == "S_plus"] <- "NVOC"
names(uk2)[names(uk2) == "day.x"] <- "day"
uk2 <- uk2[uk2$day>=446 & uk2$day <= 530, ] # start from the lowest point of S+ cases (when alpha as reached almost 100%)
#View(uk2)

################################################################################# COMPUTE RH AND RV #################################################################################

get_r <- function(tab, col_cases, col_N, col_NVOC, col_stratify = NULL){
  
  # get rH and rV from cases data, number of samples for VOC N, number of VOC NVOC
  # and if necessary stratify by some factor (typically, region)
  stopifnot(all(c(col_cases, col_N, col_NVOC, col_stratify) %in% names(tab)))
  tab$rH <- NA
  tab$rV <- NA
  # WT and variant cases:
  tab$P_H <- unlist(tab[col_cases]) * (1 - unlist(tab[col_NVOC]) / unlist(tab[col_N]))
  tab$P_V <- unlist(tab[col_cases]) * unlist(tab[col_NVOC]) / unlist(tab[col_N])
  if(!is.null(col_stratify)){
    myfactors <- unlist(unique(tab[col_stratify]))
  } else {
    col_stratify <- "identity"
    tab$identity <- 1
    myfactors <- c(1)
  }
  for(ff in myfactors){
    idx <- which(tab[col_stratify] == ff) # sub-factor
    ndays <- length(idx)
    tmp <- tab[idx, ]
    tmp[2:ndays, "rH"] <- log(unlist(tmp[2:ndays, "P_H"] / tmp[1:(ndays-1), "P_H"]))
    tmp[2:ndays, "rV"] <- log(unlist(tmp[2:ndays, "P_V"] / tmp[1:(ndays-1), "P_V"]))
    tab[idx, "rH"] <- unlist(tmp[, "rH"])
    tab[idx, "rV"] <- unlist(tmp[, "rV"])
  }
  return(tab)
}

uk <- get_r(tab = uk, col_cases = "weekly_cases", col_N = "N", col_NVOC = "NVOC")
uk2 <- get_r(tab = uk2, col_cases = "weekly_cases", col_N = "N", col_NVOC = "NVOC")
uk_region <- get_r(tab = uk_region , col_cases = "weekly_cases", col_N = "N", col_NVOC = "NVOC")
london <- get_r(tab = london, col_cases = "weekly_cases", col_N = "N", col_NVOC = "NVOC")

# data retrieved from the csv file that can be downloaded here:
# https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea
# the adress of the csv file is here:
# https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv
# on 25th August 2021

a <- read.csv("~/Downloads/data-5.csv") # The downloaded file
head(a)
a <- a[
  with(a, order(country, year_week)),
]
View(a[which(a$country == "Italy" & a$variant == "B.1.617.2" & a$source == "TESSy"),])

voc <-  "B.1.617.2"
#voc <- "B.1.1.7"
source <- "TESSy" #; source <- "GISAID"
threshold <- 50


# selected country based on inspection of the curve on TESSy (looking sigmoidally shaped)
list_cc <- c("Austria", "Belgium", "Denmark", "France", "Germany", "Greece", "Ireland", "Italy", "Netherlands", "Norway", "Sweden")
a2 <- c()
a3 <- c()

for(cc in list_cc){
  
  idx1 <- which(a$country==cc & a$variant == voc & a$source == source)
  #idx2 <- which(a$country==cc & a$variant == voc & a$source == source)
  
  if(length(idx1) > 3){
    plot(a$percent_variant[idx1], type = "o", pch = 20, ylim = c(0, 100), main = cc)
    idx1 <- idx1[1:(length(idx1)-3)]
  
    first_over_50 <- idx1[which(a$percent_variant[idx1] > 50)[1]]
    last_below_50 <- idx1[which(a$percent_variant[idx1] < 50)]; last_below_50 <- last_below_50[length(last_below_50)]
    
    date1 <- a$year_week[last_below_50]
    date2 <- a$year_week[first_over_50]
    print(c(cc, date1, date2))
    
    a2 <- rbind(a2, a[c(last_below_50, first_over_50), ])
    a3 <- rbind(a3, a[idx1, ]) # keep all data still
  }
}
a3$week <- as.numeric(sapply(a3$year_week, function(yw) strsplit(yw, split = "-")[[1]][2]))
a3$year <- as.numeric(sapply(a3$year_week, function(yw) strsplit(yw, split = "-")[[1]][1]))
a3 <- a3[a3$year == 2021, ] # restrict to 2021
write.csv(x = a3, file = "~/ownCloud/coronavirus/variant_France/delta_frequencies_EU.csv", row.names = F)


rs <- data.frame(country = list_cc, Psmoothed0 = NA, Psmoothed1 = NA, NVOC0 = NA, Ntrue0 = NA, NVOC1 = NA, Ntrue1 = NA, rH = NA, rV = NA)
for(cc in list_cc){
  idx <- which(a2$country == cc)
  
  if(length(idx) == 2){
    
    i <- which(rs$country == cc)
    
    NW <- (1 - a2$percent_variant[idx]/100) * a2$new_cases[idx]
    NM <- a2$percent_variant[idx]/100 * a2$new_cases[idx]
    
    # fill in the r table:
    rs[i, c("Psmoothed0", "Psmoothed1")] <- a2$new_cases[idx]
    rs[i, c("Ntrue0", "Ntrue1")] <- a2$number_sequenced[idx]
    rs[i, c("NVOC0", "NVOC1")] <- a2$number_detections_variant[idx]
    rs[i, "rH"] <- (1/1) * log(NW[2]/NW[1])
    rs[i, "rV"] <- (1/1) * log(NM[2]/NM[1])
    
    
  } else {
    stop("country without two entries")
  }
}

plot(rs$rH, rs$rV, pch = 20)
abline(0,1)
lm0 <- lm(rs$rV ~ rs$rH)
text(rs$rH, rs$rV + 0.05, rs$country)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
abline(lm0$coefficients[1], lm0$coefficients[2], col = "gray", lwd = 3)

summary(lm0)$r.s
summary(lm0)$coefficients
write.csv(x = rs, file = "~/ownCloud/coronavirus/variant_France/spatial_variation_EU.csv", row.names = F)


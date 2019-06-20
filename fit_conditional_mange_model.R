library(dplyr)
library(lubridate)
library(exifr)

library(runjags)



# read in coyote data
coy <- read.csv("./data/coydata_withdate.csv", stringsAsFactors = FALSE)
coy$site <- substr(coy$surveyid, 1, 8)

coy$site <- substr(coy$surveyid, 1, 8)

# drop FA13 and SU13 for now
coy <- coy[-grep('FA13|SU13', coy$surveyid),]

# bring in master data
master <- read.csv("./data/coyote_detection_data.csv", stringsAsFactors = FALSE)

# add week to the coy data
#  this will drop images that fall outside of our sampling period
coy <- inner_join(coy, master[,c('SurveyID', 'Week', 'Date')],
                  by = c('date' = 'Date', 'surveyid' = 'SurveyID'))


tmp_master <- master
tmp_master$mrow <- 1:nrow(master)

# check to see if there are detections we need to add
test <- left_join(coy, tmp_master, by= c('surveyid' = 'SurveyID', 'date' = 'Date'))
coyote_add <- test[which(test$Coyote == 0),]

master$Coyote[coyote_add$mrow] <- 1


# remove bmt0
master <- master[-which(master$StationID == "D02-BMT0"),]
# drop duplicates
master <- master[-which(duplicated(master)==TRUE),]


master$Coyote[is.na(master$Coyote)] <- -1
master$mange[is.na(master$mange)] <- -1
# summarise down to weekly detections
master_week <- master %>% group_by(SurveyID, Week) %>% 
  summarise( Coyote = max(Coyote),
             mange = max(mange),
             station = unique(substr(SurveyID,1,8)))
master_week$mange[master_week$mange == -1] <- NA
master_week$Coyote[master_week$Coyote == -1] <- NA
# c

# need to order this by season / year. This generates the appropriate vector
#  to sort by
twenty_years_of_data <- paste(c("WI", "SP", "SU", "FA"), 
                              rep(seq(10,30,1), each = 4), 
                              sep = "" )[2:14]


# add a season column to master week, giving it levels equal to the correct season/year order.
master_week$season <- factor(substr(master_week$SurveyID, 10,13),
                             levels = twenty_years_of_data)

# order by season column
master_week <- master_week[order(master_week$station, master_week$season),]

# doing this in long format
y_det <- master_week %>% 
  group_by(SurveyID) %>% 
  summarise(J = 4 - sum(is.na(Coyote)),
            y = sum(Coyote, na.rm = TRUE),
            station = unique(station),
            season = unique(season))

y_det$y[y_det$J == 0] <- NA

y_det <- y_det[order(y_det$season, y_det$station),]

# set up arrays for analysis
week_mat <- array(master_week$Coyote, dim = c(4, 13, 122))
site_mat <- array(master_week$station, dim = c(4, 13, 122))

# remove sites that have less than 2 seasons of data
togo <- rep(0, 122)
for(i in 1:122){
  na_seas <- apply(week_mat[,,i], 2, function(x) sum(is.na(x)))
togo[i] <- ifelse(sum(na_seas==4) >11, i, 0)  
}

week_mat <- week_mat[,,-togo[togo>0] ]
site_mat <- site_mat[,,-togo[togo>0] ]
sites <- unique(as.character(site_mat))

y_det <- y_det[which(y_det$station %in% sites),]
y_det$station <- factor(y_det$station, levels = sites)

# set this up in a matrix for the detection model
#  step 1. sort by site and season
coy$season <- substr(coy$surveyid, 10,13)
coy$season <- factor(coy$season, levels = twenty_years_of_data)
coy <- coy[order(coy$season, coy$site),]

# start making numeric vectors to represent the site, season, and week

# drop data that has been excluded from occupancy analysis
coy <- coy[-which(!coy$site %in% sites),]


coy$site <- factor(coy$site, levels = sites)


coy$Week <- factor(coy$Week, levels = c('week 1', 'week 2', 'week 3', 'week 4'))
coy <- coy[order( coy$season, coy$site),]


coy$sitevec <- as.numeric(factor(coy$surveyid, levels = y_det$SurveyID))

y_det$sampvec <- as.numeric(y_det$season)
y_det$fall <- 0
y_det$fall[grep('FA', y_det$season)] <- 1
y_det$winter <- 0
y_det$winter[grep('WI', y_det$season)] <- 1

# make a vector to track each season occupancy
season_tracker <- matrix(0, ncol = 2, nrow = 13)
season_tracker[,1] <- seq(1, 1339, 103)
season_tracker[,2] <- seq(103, nrow(y_det), 103)


# CONSTRUCT COVARIATES
coy$num_sev <- as.numeric(factor(coy$severity,
                                 levels = c('None',
                                            'Mild',
                                            'Moderate',
                                            'Severe'))) - 1 

#mange detection
gamma_covs <- coy[,c('blur', 'In_color','propbodyvis')]


to_1 <- which(gamma_covs$blur>=500)
gamma_covs$blur[to_1] <- 1
gamma_covs$blur[-to_1] <- 0
gamma_covs$propbodyvis <- scale(gamma_covs$propbodyvis)
gamma_covs <- cbind(1, gamma_covs)


# bring in urbanization data
urb <- read.csv("./data/raw_covariate_data.csv")
urb <- urb[which(urb$site %in% sites),]

urb$site <- factor(urb$site, levels = sites)


urb_cov <- prcomp(urb[,-1], scale. = TRUE)


urb_mat <- data.frame( site = urb$site,
                       urb1 = as.numeric(urb_cov$x[,1]) * -1,
                       urb2 = as.numeric(urb_cov$x[,2]) * -1
                       )

# merge with y_det
y_det <- left_join(y_det, urb_mat, by = c('station' = 'site'))


temps <- data.frame(season = levels(y_det$season),
                    temp = c(    47, 70, 48, 
                             18, 42, 72, 47,
                             23, 43, 74, 45,
                             20, 39 ))
temps$temp <- as.numeric(scale(temps$temp))

# merge with y_det
y_det <- left_join(y_det, temps, by = 'season')


psi_covs <- cbind(1, y_det[,c('urb1','urb2')])

rho_covs <- cbind(1, y_det[,c('temp')])

omega_covs <- cbind(1, y_det[, c('urb1','urb2', 'temp')])

ncov_psi <- ncol(psi_covs)
ncov_rho <- ncol(rho_covs)
ncov_omega <- ncol(omega_covs)
ncov_gamma <- ncol(gamma_covs)
nseason <- max(y_det$sampvec)

z_start <- y_det$y
z_start[z_start>1] <- 1
z_start[is.na(z_start)] <- 1

z <- z_start
#z[is.na(z)] <- 1

x_guess <- rep(NA, nrow(y_det)) 
x_guess[unique(coy$sitevec[coy$Mange_signs_present == 1])] <- 1


inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = as.numeric(z),
      x = x_guess,
      psi = rnorm(ncov_psi),
      rho = rnorm(ncov_rho),
      omega = rnorm(ncov_omega),
      gamma = rnorm(ncov_gamma),
      psi_ranef = rnorm(nseason),
      rho_ranef = rnorm(nseason),
      omega_ranef = rnorm(nseason),
      tau_psi = rgamma(1,1,1),
      tau_rho = rgamma(1,1,1),
      tau_omega = rgamma(1,1,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

data_list <- list( y = y_det$y,
                   J = y_det$J,
                   q = coy$Mange_signs_present,
                   psi_cov = as.matrix(psi_covs),
                   rho_cov = as.matrix(rho_covs),
                   omega_cov = as.matrix(omega_covs),
                   gamma_cov = as.matrix(gamma_covs),
                   ncov_psi = ncov_psi,
                   ncov_rho = ncov_rho,
                   ncov_omega = ncov_omega,
                   ncov_gamma = ncov_gamma,
                   nseason = nseason,
                   st = season_tracker,
                   nsite = nrow(y_det),
                   nphoto = nrow(coy),
                   sample_vec = y_det$sampvec,
                   site_vec = coy$sitevec)


start <- Sys.time()
mout3 <- run.jags(model = "./jags_script/conditional_model.R",
                 monitor = c('psi',
                             'rho',
                             'omega',
                             'gamma',
                             'sd_psi',
                             'sd_rho',
                             'sd_omega',
                             'psi_ranef',
                             'rho_ranef',
                             'omega_ranef',
                             'n_coyote',
                             'n_mange'),
                 data = data_list,
                 n.chains = 4,
                 inits = inits,
                 adapt = 1000,
                 burnin = 20000,
                 sample = ceiling(50000/4),
                 thin = 2,
                 module = 'glm',
                 method = 'parallel')
end <- Sys.time()
end - start



m2 <- as.mcmc.list(mout)
saveRDS(mout, "./results/coyote_mcmc_images.RDS")
ans <- summary(mout)
str(ans)

mm <- as.matrix(as.mcmc.list(mout), chains = TRUE)



# calculate coyote without mange
newx <- cbind(1,seq(-3.5,3.5, 0.01))

psi <- mm[,2:3]
omega <- mm[,7:8]

psi_pred <- plogis(psi %*% t(newx))
omega_pred <- plogis(omega %*% t(newx))

c_nomange <- psi_pred * (1 - omega_pred)
c_mange <- psi_pred * omega_pred

c_nomange <- t(apply(c_nomange, 2, quantile, probs = c(0.025,0.5,0.975)))
c_mange <- t(apply(c_mange, 2, quantile, probs = c(0.025,0.5,0.975)))

plot(c_nomange[,2] ~ newx[,2], xlab = "Urbanization",
     ylab = "probability of", ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2)
lines(c_nomange[,1] ~ newx[,2], lty = 2)
lines(c_nomange[,3] ~ newx[,2], lty = 2)

lines(c_mange[,2] ~ newx[,2], col = "red", lwd = 2)
lines(c_mange[,1] ~ newx[,2], col = "red", lty = 2)
lines(c_mange[,3] ~ newx[,2], col = "red", lty = 2)

legend('topright', legend = c('Coyote without mange',
                              'Coyote with mange'),
       col = c('black', 'red'), lwd = 2)

# omega pred
omega_2 <- t(apply(omega_pred, 2, quantile, probs = c(0.025,0.5,0.975)))

plot(omega_2[,2] ~ newx[,2], xlab = "Urbanization",
     ylab = "probability of a coyote having mange", 
     ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2)
lines(omega_2[,1] ~ newx[,2], lty = 2)
lines(omega_2[,3] ~ newx[,2], lty = 2)

esties <- ans[1:28,1:3]

pr_coy <- ans[1:6,2]

omega_rans <- mm[,c(43:55)] 
for(i in 1:nrow(omega_rans)){
  omega_rans[i,] <- mm[i,7] + omega_rans[i,]
}

omega_prob <- apply(plogis(omega_rans), 2, quantile,
                    probs = c(0.025,0.5,0.975))
plot(omega_prob[2,] ~ c(1:13), bty = 'l', ylim = c(0,0.75),
     pch = 20, cex = 2, ylab= 'probability of mange')

for(i in 1:13){
  lines(x = rep(i, 2),
        y = c(omega_prob[1,i], omega_prob[3,i]))
}


oc <- plogis(mm[,2])

ocm <- plogis(mm[,8])

coyote_occ <- oc * (1 - ocm)
coyote_mocc <- oc *  ocm

both <- cbind( 1 - oc, coyote_occ, coyote_mocc)

apply(both, 2, quantile, probs = c(0.025,0.5,0.975))

obs <- apply(z, 2, table)

mean(obs[2, ] / obs[1,])

b1 <- plogis(-0.43) # detection coyote 
b2 <- plogis(-1.60) # detection mange



coy_det <- c(1 - b1, b1 * (1 - b2), b1 * b2)

just_coy <- c(1 - b1, b1)

just_coy <- apply(z, 2, table)

mean(just_coy[3,] / just_coy[2,])

c1 <- plogis(-0.1 - 0.64 * seq(-2,2,0.1))
c2 <- plogis(0.16 - 0.11 * seq(-2,2,0.1))

coy_mange <- cbind(1 - c1, c1 * (1 - c2), c1 * c2)

urb <- seq(-2, 2, by = 0.1)
plot(coy_mange[,2] ~ urb, type = 'l', ylim = c(0,0.5))
lines(coy_mange[,3] ~ urb, col = "purple", lwd = 2)


library(dplyr)
library(lubridate)
library(runjags)
library(coda)

# read in the coyote image data
coy <- read.csv(
  "./data/coydata_merged_sites.csv",
  stringsAsFactors = FALSE
)

# bring in main data (i.e., the detection history data).
main <- read.csv(
  "./data/coyote_detection_data.csv",
  stringsAsFactors = FALSE
)

# add week to the coy data
#  this will drop images that fall outside of our sampling period
coy <- inner_join(
  coy,
  main[,c('SurveyID', 'Week', 'Date')],
  by = c(
    'date' = 'Date',
    'surveyid' = 'SurveyID'
  )
)

# convert down to weekly
main$Coyote[is.na(main$Coyote)] <- -1
#main$mange[is.na(main$mange)] <- -1
# summarise down to weekly detections
main_week <- main %>% group_by(SurveyID, Week) %>% 
  summarise( Coyote = max(Coyote),
             # mange = max(mange),
             station = unique(substr(SurveyID,1,8)))
#main_week$mange[main_week$mange == -1] <- NA
main_week$Coyote[main_week$Coyote == -1] <- NA
 
# remove main now that we have weekly detections
rm(main)

# need to order this by season / year. This generates the appropriate vector
#  to sort by
twenty_seasons_of_data <- paste(c("WI", "SP", "SU", "FA"), 
                              rep(seq(10,30,1), each = 4), 
                              sep = "" )[2:17]


# add a season column to main week, giving it levels equal to the correct season/year order.
main_week$season <- factor(substr(main_week$SurveyID, 10,13),
                             levels = twenty_seasons_of_data)

# order by season column
main_week <- main_week[order(main_week$station, main_week$season),]

# doing this in long format
y_det <- main_week %>% 
  group_by(SurveyID) %>% 
  summarise(J = 4 - sum(is.na(Coyote)),
            y = sum(Coyote, na.rm = TRUE),
            station = unique(station),
            season = unique(season))

y_det$y[y_det$J == 0] <- NA

y_det <- y_det[order(y_det$season, y_det$station),]

# set up arrays for analysis
week_mat <- array(main_week$Coyote, dim = c(4, 16, 122))
site_mat <- array(main_week$station, dim = c(4, 16, 122))

# remove sites that have less than 2 seasons of data
togo <- rep(0, 122)
for(i in 1:122){
  na_seas <- apply(week_mat[,,i], 2, function(x) sum(is.na(x)))
  togo[i] <- ifelse(sum(na_seas==4) >14, i, 0)  
}

week_mat <- week_mat[,,-togo[togo>0] ]
site_mat <- site_mat[,,-togo[togo>0] ]

rm(na_seas)
rm(togo)
sites <- unique(as.character(site_mat))

y_det <- y_det[which(y_det$station %in% sites),]
y_det$station <- factor(y_det$station, levels = sites)

# set this up in a matrix for the detection model
#  step 1. sort by site and season
coy$season <- substr(coy$surveyid, 10,16)
coy$season <- factor(coy$season, levels = twenty_seasons_of_data)
coy <- coy[order(coy$season, coy$site),]

# start making numeric vectors to represent the site, season, and week

# drop data that has been excluded from occupancy analysis
if(any(!coy$site %in% sites)){
  coy <- coy[-which(!coy$site %in% sites),]
}

coy$site <- factor(coy$site, levels = sites)


coy$Week <- factor(coy$Week, levels = c('week 1', 'week 2', 'week 3', 'week 4'))
coy <- coy[order( coy$season, coy$site),]


coy$sitevec <- as.numeric(factor(coy$surveyid, levels = y_det$SurveyID))

y_det$sampvec <- as.numeric(y_det$season)

# make a vector to track each season occupancy
season_tracker <- matrix(0, ncol = 2, nrow = length(twenty_seasons_of_data))
season_tracker[,1] <- seq(1, nrow(season_tracker) * length(sites), length(sites))
season_tracker[,2] <- seq(length(sites), nrow(y_det), length(sites))

rm(twenty_seasons_of_data)
# CONSTRUCT COVARIATES

#mange detection
gamma_covs <- coy[,c('blur', 'In_color','propbodyvis')]

# changing blur to binary (greater or less than mean).
to_1 <- which(gamma_covs$blur>=mean(gamma_covs$blur))
gamma_covs$blur[to_1] <- 1
gamma_covs$blur[-to_1] <- 0
gamma_covs$propbodyvis <- scale(gamma_covs$propbodyvis)
gamma_covs <- cbind(1, gamma_covs)

rm(to_1)

# bring in urbanization data

urb <- read.csv(
  "./data/model_covariates.csv",
  stringsAsFactors = FALSE
)
urb <- urb[order(urb$site),]
urb <- urb[which(urb$site %in% sites),]

urb$site <- factor(urb$site, levels = sites)


urb_cov <- prcomp(urb[,-1], scale. = TRUE)

urb_mat <- data.frame( site = urb$site,
                       urb1 = as.numeric(urb_cov$x[,1]) * -1,
                       urb2 = as.numeric(urb_cov$x[,2]),
                       stringsAsFactors = FALSE
)


urb_mat$site <- factor(urb_mat$site, levels = sites)

urb_mat <- urb_mat[order(urb_mat$site),]

# merge with y_det
y_det <- left_join(y_det, urb_mat, by = c('station' = 'site'))


temps <- data.frame(season = levels(y_det$season),
                    temp = c(    47, 70, 48, 
                                 18, 42, 72, 47,
                                 23, 43, 74, 45,
                                 20, 39, 73, 53, 
                                 16))
temps$temp <- as.numeric(scale(temps$temp))

# merge with y_det
y_det <- suppressWarnings(left_join(y_det, temps, by = 'season'))


psi_covs <- cbind(1, y_det[,c('urb1','urb2')])
psi_covs <- cbind(psi_covs, psi_covs[,2] * psi_covs[,3])

rho_covs <- cbind(1, y_det[,c('temp')])

omega_covs <- cbind(1, y_det[, c('urb1','urb2')])#, 'temp')])
omega_covs <- cbind(omega_covs, omega_covs[,2]  * omega_covs[,3])
#omega_covs[which(y_det$fall == 1),4] <- 1

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


inits_resid <- function(chain){
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
      psi_resid_tau = rgamma(1,1,1),
      rho_resid_tau = rgamma(1,1,1),
      ome_resid_tau = rgamma(1,1,1),
      gam_resid_tau = rgamma(1,1,1),
      psi_resid = rnorm(data_list$nsite, 0, 0.1),
      rho_resid = rnorm(data_list$nsite, 0, 0.1),
      ome_resid = rnorm(data_list$nsite, 0, 0.1),
      gam_resid = rnorm(data_list$nphoto, 0, 0.1),
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
inits_resid2 <- function(chain){
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
      psi_resid_tau = rgamma(1,1,1),
      rho_resid_tau = rgamma(1,1,1),
      ome_resid_tau = rgamma(1,1,1),
      gam_resid_tau = rgamma(1,1,1),
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


rm(sites)
rm(site_mat)
rm(week_mat)
rm(i)
rm(z_start)
rm(temps)
rm(gamma_covs)
rm(omega_covs)
rm(psi_covs)
rm(rho_covs)
#rm(urb)

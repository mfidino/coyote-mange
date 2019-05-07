library(dplyr)
library(lubridate)
library(exifr)

# read in coyote data
coy <- read.csv("./data/coydata_withdate.csv", stringsAsFactors = FALSE)


# get whether or not mange was detected during each surveyID
mange_per_survey <- coy %>% 
  group_by(surveyid) %>% 
  summarise(mange = max(Mange_signs_present))

# bring in master data
master <- read.csv("./data/coyote_detection_data.csv", stringsAsFactors = FALSE)

# add mange detections to the master
master <- left_join(master, mange_per_survey, by = c("SurveyID" = "surveyid"))

# remove bmt0
master <- master[-which(master$StationID == "D02-BMT0"),]
# drop duplicates
master <- master[-which(duplicated(master)==TRUE),]

# make a column to indicate the state of the detection
master$state <- NA
# for any time when camera is operable put state 1
#  this indicates that coyote were not detected
master$state[!is.na(master$Coyote)] <- 1
# for any time when a coyote was detected put state 2
#  this indicates that a coyote was detected
master$state[which(master$Coyote == 1)] <- 2
# for any time when mange was detected put state 3
#  this indicates that a coyote with mange was detected
master$state[which(master$mange == 1)] <- 3
# return NA values when Coyote = NA
master$state[which(is.na(master$Coyote))] <- NA

# start condensing down to weekly detections
#  step 1: put a -1 in this column if camera not operational
master$weekly <- is.na(master$state) * 0 + -1
#  step 2: put the state column in when there is data
master$weekly[-which(is.na(master$state))] <- master$state[-which(is.na(master$state))]

# summarise down to weekly detections
master_week <- master %>% group_by(SurveyID, Week) %>% 
  summarise(state = max(weekly), station = unique(StationID))

# need to order this by season / year. This generates the appropriate vector
#  to sort by
twenty_years_of_data <- paste(c("WI", "SP", "SU", "FA"), 
                              rep(seq(10,30,1), each = 4), 
                              sep = "" )
# add a season column to master week, giving it levels equal to the correct season/year order.
master_week$season <- factor(substr(master_week$SurveyID, 10,13), levels = twenty_years_of_data)

# order by season column
master_week <- master_week[order(master_week$station, master_week$season),]

# put an NA in place when state is -1
master_week$state[master_week$state == -1] <- NA

# set up arrays for analysis
week_mat <- array(master_week$state, dim = c(4, 13, 122))
site_mat <- array(master_week$station, dim = c(4, 13, 122))

# remove sites that have less than 2 seasons of data
togo <- rep(0, 122)
for(i in 1:122){
  na_seas <- apply(week_mat[,,i], 2, function(x) sum(is.na(x)))
togo[i] <- ifelse(sum(na_seas==4) >11, i, 0)  
}

week_mat <- week_mat[,,-togo[togo>0] ]
sites <- unique(as.character(site_mat[,,which(togo==0)]))

# get maximum state for z initial values
z_start <- apply(week_mat, c(3,2), max, na.rm = TRUE)
z_start[is.infinite(z_start)] <- NA

# set this up in a matrix for the detection model
#  step 1. sort by site and season
coy$season <- substr(coy$surveyid, 10,13)
coy$season <- factor(coy$season, levels = twenty_years_of_data)
coy <- coy[order(coy$season, coy$site),]


detect_covs <- left_join(longshot, prop_clear, by = "SurveyID")
detect_covs[is.na(detect_covs)] <- 0

detect_covs <- detect_covs[-which(duplicated(detect_covs$SurveyID)==TRUE),]
mu_mat <- t(matrix(detect_covs$mu, nrow = 13, ncol =  122)[,-togo[togo>0]])
sigma_mat <- t(matrix(detect_covs$sigma, nrow = 13, ncol =  122)[,-togo[togo>0]])
inx_mat <- t(matrix(detect_covs$inx, nrow = 13, ncol =  122)[,-togo[togo>0]])

# bring in n photos of canids

canid_photo <- read.csv("n_photo_canids.csv", stringsAsFactors = FALSE)

canid_photo <- canid_photo %>% group_by(SurveyID) %>% 
  summarise(canids = sum(mammals, na.rm = TRUE))
scaled_canid <- scale(canid_photo$canids)
canid_photo$canids <- scaled_canid
scale_covs$scaled_canids <- scaled_canid

canid <- left_join(longshot, canid_photo, by = "SurveyID")

canid <- canid[-which(duplicated(canid$SurveyID) == TRUE),]
canid_mat <- t(matrix(canid$canids, nrow = 13, ncol = 122)[,-togo[togo>0]])
canid_mat[is.na(canid_mat)] <- 0

# bring in urbanization data
urb <- read.csv("raw_covariate_data.csv")
urb <- urb[which(urb$site %in% sites),]

urb <- prcomp(urb[,-1], scale. = TRUE)


urb_mat <- as.matrix(urb$x[,1]) * -1

library(runjags)
library(rjags)

week_mat <- aperm(week_mat, c(3,2,1))


acovs <- matrix(urb_mat)
bcovs <- array(NA, dim = c(102,13,2))

bcovs[,,1] <- array(urb_mat, dim = c(102, 13))
bcovs[,,2] <- canid_mat

dcovs <- acovs

fcovs <- array(NA, dim = c(102,13,4))

fcovs[,,1] <- mu_mat
fcovs[,,2] <- sigma_mat
fcovs[,,3] <- inx_mat
fcovs[,,4] <- urb_mat

sea_vec <- as.character(unique(longshot$season))
sea_vec <- factor(substr(sea_vec, 1, 2), levels = c("SP", "SU", "FA", "WI"))
sea_vec <- as.numeric(sea_vec)

year_vec <- as.character(unique(longshot$season))                  
year_vec <- factor(substr(year_vec, 3, 4))
year_vec <- as.numeric(year_vec)

# make inits for model

z <- matrix(NA, nrow = 102, ncol = 13)

for(i in 1:102){
  for(j in 1:13){
  z[i,j] <- max(week_mat[i,j,], na.rm = TRUE)
  }
}
z[is.infinite(z)]  <- 0
z <- z+1
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = as.matrix(z),
      a0 = rnorm(1),
      a = rnorm(1),
      aseasonpar = rnorm(3),
      b0 = rnorm(1),
      b = rnorm(2),
      bseasonpar = rnorm(3),
      d0 = rnorm(1),
      d = rnorm(1),
      dseasonpar = rnorm(3),
      f0 = rnorm(1),
      f = rnorm(4),
      fseasonpar = rnorm(3),
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

week_mat[!is.na(week_mat)] <- week_mat[!is.na(week_mat)] + 1

data_list <- list( y = week_mat, acovs = acovs, bcovs = bcovs, dcovs = dcovs, 
           fcovs = fcovs, sea_vec = sea_vec,
           nyear = 13, nsite = 102, nsurvey = 4, nstate = 3)

library(rjags)
library(runjags)

load.module("glm")

mout <- run.jags(model = "./jags_script/new_coyote_mange_model.R",
                 monitor = c("a0", "a", "aseason",
                             "b0", "b", "bseason",
                             "d0", "d", "dseason",
                             "f0", "f", "fseason"),
                 data = data_list,
                 n.chains = 7,
                 inits = inits,
                 adapt = 1000,
                 burnin = 20000,
                 sample = ceiling(25000/7),
                 thin = 3,
                 summarise = FALSE,
                 plots = FALSE,
                 method = "parallel")
saveRDS(mout, "./results/coyote_mcmc.RDS")
ans <- summary(mout)
str(ans)

mm <- as.matrix(as.mcmc.list(mout), chains = TRUE)

esties <- ans[1:28,1:3]

pr_coy <- ans[1:6,2]

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


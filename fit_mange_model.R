library(dplyr)
library(lubridate)
library(exifr)

# read in coyote data
coy <- read.csv("./data/coydata_withdate.csv", stringsAsFactors = FALSE)
coy$site <- substr(coy$surveyid, 1, 8)

#season_swapper <- function(season){
#  switch(season,
#         'Spring' = "SP",
#         'Winter' = 'WI',
#         'Fall' = 'FA',
#         'Summer' = 'SU')
#}
#
## generate survey ID
#coy$surveyid <- paste0(substr(coy$new_file_name, 1, 8),"-",
#                       unname(sapply(coy$Season, season_swapper)),
#                       substr(coy$Year, 3,4))
coy$site <- substr(coy$surveyid, 1, 8)

# drop FA13 and SU13 for now
coy <- coy[-grep('FA13|SU13', coy$surveyid),]

#to_combine <- read.table("./data/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")
#
#o1 <- to_combine[to_combine$group==1,]
#o1$site_no_number <- paste0(o1$site_no_number, "0")
#
#for(i in 1:nrow(o1)){
#  if(any(coy$site == o1$Site[i])){
#  coy$surveyid[coy$site == o1$Site[i]] <- 
#    paste0(o1$site_no_number[i], "-", substr(coy$surveyid[coy$site == o1$Site[i]], 10,13))
#  }
#}
#coy$site <- substr(coy$surveyid, 1, 8)

# get the date / time from each image
#  step 1: construct the file path
#tmp_path <- paste0("X:/PEOPLE/Mason Fidino/coyote photos/output_",
#                   substr(coy$surveyid, 10, 13), "/", coy$new_file_name)
#
#test <- read_exif(tmp_path, tags =c('DateTimeOriginal'))
#
#coy$date <- sapply(strsplit(test$DateTimeOriginal," "), function(x) x[1]) %>% 
#  gsub(':',"-", .)

# convert D05-MCH to D05-MHC
#coy$surveyid <- gsub('D05-MCH', 'D05-MHC', coy$surveyid)


# check which coy not in master



#write.csv(coy, './data/coydata_withdate.csv', row.names = FALSE)




# bring in master data
master <- read.csv("./data/coyote_detection_data.csv", stringsAsFactors = FALSE)



# add week to the coy data
#  this will drop images that fall outside of our sampling period
coy <- inner_join(coy, master[,c('SurveyID', 'Week', 'Date')],
                  by = c('date' = 'Date', 'surveyid' = 'SurveyID'))

# get whether or not mange was detected during each surveyID
mange_per_survey <- coy %>% 
  group_by(surveyid, Week) %>% 
  summarise(mange = max(Mange_signs_present))



#to_check <- which(!coy$surveyid %in% master$SurveyID)
#to_check <- unique(coy$surveyid[to_check])
#my_surveys <- unique(master$SurveyID)
#to_check
#unique(my_surveys[grep('D02', my_surveys)])
#
#coy$surveyid <- gsub('S03-LCP', 'S03-LPC', coy$surveyid)
#master$SurveyID <- gsub('S03-LCP', 'S03-LPC', master$SurveyID)
#master$SurveyID <- gsub('D05-MCH', 'D05-MHC', master$SurveyID)

# add mange detections to the master
master <- left_join(master, mange_per_survey, by = c("SurveyID" = "surveyid",
                                                     'Week'))

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
  summarise(state = max(weekly), station = unique(substr(SurveyID,1,8)))

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
surv_mat <- array(master_week$SurveyID, dim = dim(site_mat))

# remove sites that have less than 2 seasons of data
togo <- rep(0, 122)
for(i in 1:122){
  na_seas <- apply(week_mat[,,i], 2, function(x) sum(is.na(x)))
togo[i] <- ifelse(sum(na_seas==4) >11, i, 0)  
}

week_mat <- week_mat[,,-togo[togo>0] ]
site_mat <- site_mat[,,-togo[togo>0] ]
sites <- unique(as.character(site_mat))

# get maximum state for z initial values
z_start <- apply(week_mat, c(3,2), max, na.rm = TRUE)
z_start[is.infinite(z_start)] <- 1

# set this up in a matrix for the detection model
#  step 1. sort by site and season
coy$season <- substr(coy$surveyid, 10,13)
coy$season <- factor(coy$season, levels = twenty_years_of_data)
coy <- coy[order(coy$season, coy$site),]

# start making numeric vectors to represent the site, season, and week

# drop data that has been excluded from occupancy analysis
coy <- coy[-which(!coy$site %in% sites),]
coydet <- coy


coy$site <- factor(coy$site, levels = sites)
coy$Week <- factor(coy$Week, levels = c('week 1', 'week 2', 'week 3', 'week 4'))
coy <- coy[order(coy$site, coy$season, coy$Week),]

coydet <- coy





coy_summary <- coydet %>% group_by(site, season, Week) %>% 
  summarise(n_photo = length(new_file_name))

coy_summary$has_photo <- 1

# the vectors that indicate where the photos are:
coy_summary$min <- 0
coy_summary$max <- 0
coy_summary$min[1] <- 1
coy_summary$max[1] <- 11
for(i in 2:nrow(coy_summary)){
  coy_summary$min[i] <- coy_summary$max[i-1] +1
  coy_summary$max[i] <- coy_summary$min[i] + (coy_summary$n_photo[i] - 1)
}
colnames(coy_summary) <- tolower(colnames(coy_summary))

test <- expand.grid(site = sites, season = unique(master_week$season),
                    week = unique(master_week$Week))
test <- test[order(test$site, test$season, test$week),]

hm <- left_join(test, coy_summary, by = c('site', 'season', 'week'))
hm$has_photo[is.na(hm$has_photo)] <- 0
hm$min[is.na(hm$min)] <- 1
hm$max[is.na(hm$max)] <- 1

hm <- hm[order( hm$site,hm$season),]


coydet$site <- as.numeric(coydet$site)
coydet$season <- as.numeric(coydet$season)
coydet$Week <- as.numeric(coydet$Week)

hp <- array(hm$has_photo, dim = dim(site_mat))
hp[is.na(hp)] <- 0
mmin <- array(hm$min, dim = dim(site_mat))
mmin[is.na(mmin)] <- 1
mmax <- array(hm$max, dim = dim(site_mat))
mmax[is.na(mmax)] <- 1

mange_covs <- coy[,c('blur', 'In_color', 'Season','propbodyvis')]
mange_covs$Season[mange_covs$Season!='Winter'] <- 0
mange_covs$Season[mange_covs$Season=='Winter'] <- 1
mange_covs$Season <- as.numeric(mange_covs$Season)

to_1 <- which(mange_covs$blur>=500)
mange_covs$blur[to_1] <- 1
mange_covs$blur[-to_1] <- 0
mange_covs$propbodyvis <- scale(mange_covs$propbodyvis)



# toy_list <- list(mange_signs_present = coy$Mange_signs_present,
#                  mange_covs = as.matrix(mange_covs),
#                  mmin = mmin,
#                  mmax = mmax,
#                  hp = hp,
#                  nsite = 103,
#                  nyear = 13,
#                  nweek = 4,
#                  I = ncol(mange_covs),
#                  nphoto = nrow(coy))
# 
# library(runjags)
# 
# mout <- run.jags(model = "./jags_script/test_binomial.R",
#                  monitor = c(
#                              "f0", "f"),#,'mange_mu', 'pr_mangy_coyote'),
#                  data = toy_list,
#                  n.chains = 1,
#                  adapt = 10000,
#                  burnin = 10000,
#                  sample = 20000,
#                  thin = 1,
#                  summarise = FALSE,
#                  plots = FALSE,
#                  method = 'parallel')
# 
# hey <- as.matrix(as.mcmc.list(mout), chains = TRUE)
# 
# summary(mout, vars = 'f')
# 
# mmu <- hey[,grep('mu', colnames(hey))]
# mco <- hey[,grep('coyote', colnames(hey))]
# colnames(mco)[113]
# mco[,113]
# 
# detect_covs <- left_join(longshot, prop_clear, by = "SurveyID")
# detect_covs[is.na(detect_covs)] <- 0

#detect_covs <- detect_covs[-which(duplicated(detect_covs$SurveyID)==TRUE),]
#mu_mat <- t(matrix(detect_covs$mu, nrow = 13, ncol =  122)[,-togo[togo>0]])
#sigma_mat <- t(matrix(detect_covs$sigma, nrow = 13, ncol =  122)[,-togo[togo>0]])
#inx_mat <- t(matrix(detect_covs$inx, nrow = 13, ncol =  122)[,-togo[togo>0]])

# bring in n photos of canids

#canid_photo <- read.csv("n_photo_canids.csv", stringsAsFactors = FALSE)

#canid_photo <- canid_photo %>% group_by(SurveyID) %>% 
#  summarise(canids = sum(mammals, na.rm = TRUE))
#scaled_canid <- scale(canid_photo$canids)
#canid_photo$canids <- scaled_canid
#scale_covs$scaled_canids <- scaled_canid

#canid <- left_join(longshot, canid_photo, by = "SurveyID")

#canid <- canid[-which(duplicated(canid$SurveyID) == TRUE),]
#canid_mat <- t(matrix(canid$canids, nrow = 13, ncol = 122)[,-togo[togo>0]])
#canid_mat[is.na(canid_mat)] <- 0

# bring in urbanization data
urb <- read.csv("./data/raw_covariate_data.csv")
urb <- urb[which(urb$site %in% sites),]

urb <- prcomp(urb[,-1], scale. = TRUE)


urb_mat <- as.matrix(urb$x[,1]) * -1

#library(runjags)
#library(rjags)

week_mat <- aperm(week_mat, c(3,2,1))
mmin <- aperm(mmin, c(3,2,1))
mmax <- aperm(mmax, c(3,2,1))
hp <- aperm(hp, c(3,2,1))

#mmin <- mmin[-togo[togo>0],,]
#mmax <- mmax[-togo[togo>0],,]
#hp <- hp[-togo[togo>0],,]


acovs <- bcovs <- matrix(urb_mat)

#bcovs <- array(NA, dim = c(103,13,2))

#bcovs[,,1] <- array(urb_mat, dim = c(103, 13))
#bcovs[,,2] <- canid_mat

dcovs <- acovs

#fcovs <- array(NA, dim = c(102,13,4))

#fcovs[,,1] <- mu_mat
#fcovs[,,2] <- sigma_mat
#fcovs[,,3] <- inx_mat
#fcovs[,,4] <- urb_mat

sea_vec <- substr(as.character(unique(master_week$season)), 1,2)
sea_vec <- factor(sea_vec, levels = c("SP", "SU", "FA", "WI"))
sea_vec <- as.numeric(sea_vec)

#year_vec <- as.character(unique(longshot$season))                  
#year_vec <- factor(substr(year_vec, 3, 4))
#year_vec <- as.numeric(year_vec)

# make inits for model

#z <- matrix(NA, nrow = 102, ncol = 13)

#for(i in 1:102){
#  for(j in 1:13){
#  z[i,j] <- max(week_mat[i,j,], na.rm = TRUE)
#  }
#}
#z[is.infinite(z)]  <- 0
#z <- z+1
z <- z_start
#z[is.na(z)] <- 1
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = as.matrix(z),
      a0 = rnorm(1,-2),
      a = rnorm(1,-2),
      aseasonpar = rnorm(3,-2),
      b0 = rnorm(1,-2),
      b = rnorm(1,-2),
      bseasonpar = rnorm(3,-2),
      d0 = rnorm(1,-2),
      d = rnorm(1,-2),
      dseasonpar = rnorm(3,-2),
      f0 = rnorm(1,-5),
      f = rnorm(4,-1),
      fseasonpar = rnorm(3,-2),
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

#week_mat[!is.na(week_mat)] <- week_mat[!is.na(week_mat)] + 1


#toy_list <- list(mange_signs_present = coy$Mange_signs_present,
#                 mange_covs = as.matrix(mange_covs),
#                 mmin = mmin,
#                 mmax = mmax,
#                 hp = hp,
#                 nsite = 103,
#                 nyear = 13,
#                 nweek = 4,
#                 I = ncol(mange_covs),
#                 nphoto = nrow(coy))

data_list <- list( y = week_mat, acovs = acovs, bcovs = bcovs, dcovs = dcovs, 
                   gcovs = acovs,
            sea_vec = sea_vec,
           nyear = 13, nsite = 103, nsurvey = 4, nstate = 3,
           mange_signs_present = coy$Mange_signs_present,
           mange_covs = as.matrix(mange_covs),
           mmin = mmin,
           mmax = mmax,
           hp = hp,
           nphoto = nrow(coy),
           nweek=4)

library(rjags)
library(runjags)

load.module("glm")

mout <- run.jags(model = "./jags_script/coyote_mange_model_with_images_intercept.R",
                 monitor = c("a0", "a", "aseason",
                             "b0", "b", "bseason",
                             "d0", "d", "dseason",
                             "f0", "f", "fseason",
                             'g0', 'g'),
                 data = data_list,
                 n.chains = 7,
                 inits = inits,
                 adapt = 1000,
                 burnin = 40000,
                 sample = ceiling(100000/7),
                 thin = 3,
                 summarise = FALSE,
                 plots = FALSE,
                 method = "parallel")

m2 <- as.mcmc.list(mout)
saveRDS(mout, "./results/coyote_mcmc_images.RDS")
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


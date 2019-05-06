library(dplyr)
library(lubridate)

coy <- read.csv("coyote_mange_data.csv", stringsAsFactors = FALSE)

coy$Mange_signs_present[is.na(coy$Mange_signs_present)] <- 0
# remove no coyote
#coy <- coy[-which(coy$Not_coyote==1),]

colnames(coy)  <- gsub("\\.", "_", colnames(coy))

colnames(coy) <- tolower(colnames(coy))


# remove some columns

togo <- c(1,2,  4, 7, 9, 11:14, 16:20)

coy <- coy[,-togo]
csv_files <- list.files("T:/PEOPLE/Mason Fidino/coyote photos/", ".csv",
                        full.names = TRUE)

wi11 <- list.files("T:/PEOPLE/Mason Fidino/coyote photos/output_WI11", ".JPG",
                   full.names = TRUE)
just_names <- list.files("T:/PEOPLE/Mason Fidino/coyote photos/output_WI11", ".JPG")

hm <- exifr(wi11, exiftoolargs = "-DateTimeOriginal")
hm$SourceFile <- just_names

hm$DateTimeOriginal <- substr(hm$DateTimeOriginal, 1, 10)
hm$DateTimeOriginal <- gsub(":", "-",hm$DateTimeOriginal )
hm$DateTimeOriginal <- ymd(hm$DateTimeOriginal)

colnames(hm) <- c("new_file_name", "date")




season_list <- vector("list", length = length(csv_files))

for(i in 1:length(csv_files)){
  my_csv <- read.csv(csv_files[i], stringsAsFactors = FALSE)
  my_csv$date <- mdy(my_csv$date)  
  season_list[[i]] <- left_join(my_csv[,c(4,6)], coy, by = "new_file_name" )
}

season_list[[15]] <- left_join(hm, coy, by = "new_file_name")
test <- bind_rows(season_list)

coy <- coy[-which(coy$not_coyote == 1),]
test <- test[-which(test$not_coyote == 1),]
test <- test[-which(is.na(test$year)),]

coy <- test

coy$site <- substr(coy$new_file_name, 1, 8)

coy$date <- ymd(coy$date)

new_files <- list.files("T:/PEOPLE/Mason Fidino/coyote photos/mystery_photos/Maureen's edited files/Mystery photos/mystery coyotes", ".JPG|.jpg", full.names = TRUE)

library(exifr)

file_times <- exifr(new_files, exiftoolargs = "-DateTimeOriginal")

file_times$Date <- gsub(":", "/", substr(file_times$DateTimeOriginal, 1 ,10))
file_times$Date <- ymd(file_times$Date)

file_times$Time <- substr(file_times$DateTimeOriginal, 12, 20)

file_sites <- strsplit(file_times$SourceFile, "/")
file_sites <- sapply(file_sites, function(x) x[length(x)])
file_times$site <- substr(file_sites, 1, 8)

dt <- paste(coy$date, coy$time, sep = "-")
dt2 <- paste(file_times$Date, file_times$Time, sep = "-")

file_to_remove <- which(dt2 %in% dt)

file_times <- file_times[-file_to_remove,]

new_coyote_photos <- read.csv("new_coyote_photos.csv", stringsAsFactors = FALSE)

# change the source file stuff a bit so it matches
new_coyote_photos$path <- gsub("\\\\", "/", new_coyote_photos$path )
new_coyote_photos$path <- gsub("P:", "T:",new_coyote_photos$path )

site_stuff <- strsplit(new_coyote_photos$path, "/")
site_stuff <- sapply(site_stuff, function(x) x[length(x)])
new_coyote_photos$Site <- substr(site_stuff, 1, 8)

# only keep those in file times
new_coyote_photos <- new_coyote_photos[which(new_coyote_photos$path %in% 
                                               file_times$SourceFile), ]
new_coyote_photos <- new_coyote_photos[order(new_coyote_photos$path),]
file_times <- file_times[order(file_times$SourceFile),]


new_coyote_photos$date <- file_times$Date
new_coyote_photos$time <- file_times$Time
# remove not coyote
new_coyote_photos <- new_coyote_photos[- which(new_coyote_photos$Not_coyote == 1),]

colnames(new_coyote_photos)  <- gsub("\\.", "_", colnames(new_coyote_photos))

colnames(new_coyote_photos) <- tolower(colnames(new_coyote_photos))


colnames(new_coyote_photos)

new_coyote_photos <- new_coyote_photos[,-c(1, 3, 5, 7, 8, 9:10,
                                           12, 13, 14, 15:16)]
coy <- coy[,-which(colnames(coy)=="species")]
new_coyote_photos <- new_coyote_photos[,colnames(coy)]
new_coyote_photos$mange_signs_present[is.na(new_coyote_photos$mange_signs_present)] <- 0

sea <- strsplit(new_coyote_photos$new_file_name, "-")
sea <- sapply(sea, function(x) x[length(x)])
surid <- substr(sea, 1, 4)

new_coyote_photos$SurveyID <- paste(new_coyote_photos$site, surid, sep = "-")

coy$SurveyID <- paste0(substr(coy$new_file_name,1, 9 ), 
                       toupper(substr(coy$season, 1, 2)),
                       substr(coy$year, 3, 4))

p1 <- paste(coy$site, coy$date, coy$time)
p2 <- paste(new_coyote_photos$site,new_coyote_photos$date, new_coyote_photos$time )

which(p1 %in% p2)
which(p2 %in% p1)
coy$site <- substr(coy$SurveyID, 1, 8)



coy <- coy[-which(duplicated(coy)== TRUE),]

lol <- rbind(coy, new_coyote_photos)

p1 <- paste(lol$site, lol$date, lol$time)

sum(duplicated(p1))

lol <- lol[-which(duplicated(p1)==TRUE),]

coy <- lol
colnames(coy)[4] <- "state"




to_combine <- read.table("./../Transect/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

for(i in 1:nrow(o1)){
  if(any(coy$site == o1$Site[i])){
  coy$SurveyID[coy$site == o1$Site[i]] <- 
    paste0(o1$site_no_number[i], "-", substr(coy$SurveyID[coy$site == o1$Site[i]], 10,13))
  }
}
coy$site <- substr(coy$SurveyID, 1, 8)

coy_daily <- coy %>% group_by(SurveyID, date) %>% 
  summarise(state = as.numeric(sum(state)>0))

# bring in master data

master <- read.csv("coyote_detection_data.csv", stringsAsFactors = FALSE)

master$joiny <- paste(master$SurveyID, master$Date, sep = "-")
coy_daily$joiny <- paste(coy_daily$SurveyID, coy_daily$date, sep = "-")

sid <- substr(coy$SurveyID, 10,13)

cheeky <- which(sid %in% c("SU13", "FA13"))

lj <- left_join(master, coy_daily[,c(3:4)], by = "joiny")

# remove bmt0
lj <- lj[-which(lj$StationID == "D02-BMT0"),]

lj <- lj[-which(duplicated(lj)==TRUE),]

hm <- rle(lj$SurveyID)

table(lj$state[is.na(lj$Coyote)])

lj$Coyote[!is.na(lj$state)] <- 1
lj$state[lj$Coyote== 1 & is.na(lj$state)]

to_1 <- which(lj$Coyote == 1 & is.na(lj$state) == TRUE)

#lj$state[to_1] <- 0


lj$state[lj$state==1] <- 2
lj$state[lj$state==0] <- 1
lj$state[lj$Coyote == 0] <- 0

lj$weekly <- is.na(lj$state) * 0 + -1
lj$weekly[-is.na(lj$state)] <- lj$state[-is.na(lj$state)]

longshot <- lj %>% group_by(SurveyID, Week) %>% 
  summarise(state = max(weekly), station = unique(StationID))

twenty_years_of_data <- paste(c("WI", "SP", "SU", "FA"), 
                              rep(seq(10,30,1), each = 4), 
                              sep = "" )
longshot$season <- factor(substr(longshot$SurveyID, 10,13), levels = twenty_years_of_data)

  longshot <- longshot[order(longshot$station, longshot$season),]

week_mat <- array(longshot$state, dim = c(4, 13, 122))
site_mat <- array(longshot$station, dim = c(4, 13, 122))

# remove sites that have less than 2 seasons of data
togo <- rep(0, 122)
for(i in 1:122){
  na_seas <- apply(week_mat[,,i], 2, function(x) sum(is.na(x)))
togo[i] <- ifelse(sum(na_seas==4) >11, i, 0)  
}

week_mat <- week_mat[,,-togo[togo>0] ]
sites <- unique(as.character(site_mat[,,which(togo==0)]))

# make the detection covariates

prop_clear <- coy[coy$site %in% sites,]

prop_clear <- prop_clear %>% group_by(SurveyID) %>% 
  summarise(a = sum(clear_photo)+1,
            b = length(clear_photo) - sum(clear_photo) + 1) %>% 
  group_by(SurveyID) %>% 
  summarise(mu = a / (a + b),
            sigma = (a * b) / ((a + b)^2 * (a + b + 1)))
prop_clear$inx <- prop_clear$mu * prop_clear$sigma

scale_covs <-list()
scaled_prop <- scale(prop_clear[,c("mu", "sigma", "inx")])
prop_clear[,c("mu", "sigma", "inx")] <- data.frame(scaled_prop)
scale_covs$scaled_prop <- scaled_prop

# set this up in a matrix for the detection model

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

mout <- run.jags(model = "new_coyote_mange_model.R",
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
saveRDS(mout, "coyote_mcmc.RDS")
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


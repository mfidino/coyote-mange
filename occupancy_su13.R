coy <- read.csv("./data/coydata_to_WI14.csv", stringsAsFactors = FALSE)

# fix station ids
to_combine <- read.table("./data/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

# SU13 goes from 7/2/2013 to 7/29/2013

su <- read.csv("./data/coyote_su13.csv", stringsAsFactors = FALSE)

for(i in 1:nrow(o1)){
  
  su$Site <-  gsub(paste0(o1$Site[i],'(.*)'), 
                          paste0(o1$site_no_number[i],"\\1") , 
                          su$Site)
  
}

su <- su[which(su$Site %in% master$StationID),]

su$SurveyID <- paste0(su$Site,"-SU13")
su$Season <- "SU13"

# fold in the 7 days, counting to make sure there

yoo <- function(x){
  ifelse(sum(is.na(x) == 2) | length(x)==1, sum(x), sum(x, na.rm = TRUE))
}

su <- su %>% group_by(Site, Season, SurveyID) %>% 
  summarise_if(  is.integer, yoo)

# add more sites

su <- reshape2::melt(su, id.vars = c("Site", "Season", "SurveyID"))

su <- su[order(su$Site),]

colnames(su) <- c("StationID", "Season", "SurveyID", "Date", "Coyote")

levels(su$Date) <- as.character(seq(ymd("2013-7-2"), ymd("2013-7-29"), "1 day"))
su$Date <- as.character(su$Date)

su$IDWeek <- paste(su$SurveyID, su$Date)

su$Week <- rep(rep(paste("Week", 1:4), each = 7), 123)
class(master$Date)

su$SeasonWeek <- paste(su$Season, su$Week)


su <- su[,colnames(master)]

which(colnames(su) %in% colnames(master))

write.csv(su, "./data/su13_master.csv", row.names = FALSE)

# now we've got to prepare the coyote image stuff

library(uwinutils)

sites <- SELECT(
  "SELECT * FROM CameraLocations cl WHERE cl.areaID = 1"
)

coy_photo <- read.csv("./../../../Desktop/coyote_2013/detections_to_update_2019-08-27.csv",
                      stringsAsFactors = FALSE)

coy_photo <- coy_photo[,c("locationName", "photoName")]

coy_photo <- left_join(coy_photo, sites[,c(2,6)], by = c("locationName" = "fullName"))

coy_b <- read.csv("./../../../Desktop/coyote_2013/blur_output.csv",
                  stringsAsFactors = FALSE)

coy_b$path <- strsplit(
  coy_b$path,
  "\\\\"
 ) %>% 
  sapply(
    .,
    "[[",
    4
  )


coy_photo <- left_join(coy_photo, coy_b, by = c("photoName" = "path"))
coy_photo$File_name_order <- NA
coy_photo$Mange_signs_present <- NA
coy_photo$In_color <- NA
coy_photo$Season <- "Summer"
coy_photo$Year <- 2013
coy_photo$Site <- substr(coy_photo$locationAbbr, 5,7)
coy_photo$propbodyvis <- NA
coy_photo$propbodyvismange <- NA
coy_photo$severity <- NA
coy_photo$confidence <- ""
coy_photo$surveyid <- paste0(coy_photo$locationAbbr,"-SU13")

# get the date photos taken

tmp <- list.files("../../../Desktop/coyote_2013/","*.jpg",  full.names = TRUE)
ptm <- read_exif(tmp, "DateTimeOriginal")

tmp <- strsplit(ptm$DateTimeOriginal, " ") %>% sapply(., "[[", 1)
tmp <- gsub(":", "-", tmp)

tmp2 <- list.files("../../../Desktop/coyote_2013/","*.jpg")

to_join <- tibble(photoName = tmp2, date = tmp)

coy_photo <- left_join(coy_photo, to_join, by = "photoName")

coy_photo <- coy_photo[,-1]

colnames(coy_photo)[1:3] <- c("new_file_name", "site", "blur")

coy_photo <- coy_photo[,colnames(coy)]
coy_photo <- coy_photo[-which(duplicated(coy_photo$new_file_name)),]

write.csv(coy_photo, "su13_photos_formatted.csv", row.names = FALSE)

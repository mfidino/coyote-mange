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


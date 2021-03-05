###################################
#
# Plotting mange on the landscape
# 
# Written by M. Fidino 1/8/2020 mdy
#
###################################

library(dplyr)
library(data.table)
library(sf)
library(prettymapr)


# Note: There are shapefiles that are hardcoded in here, which I got
#  from the Illinois geospatial clearinghouse
##############
# I use this a lot to get everything into the same
#  spatial coordinates, so I'm making it an object.
utm_crs <- 26916
##############

tree_imperv <- data.table::fread(
  "./data/map_data/tree_imperv_mapdata.csv",
  data.table = FALSE
)

# sort it
tree_imperv <- tree_imperv[order(tree_imperv$LocationName),]

tree_imperv$num <- as.integer(tree_imperv$LocationName)

# read in housing density
housing <- read.csv(
  "./data/map_data/coyote_mange_hu10.csv",
  stringsAsFactors = FALSE
)

# read in the gridpoints
gridpoint <- sf::st_read(
  "./data/gridpoints.shp",
  crs = utm_crs
)

# convert name to simplify join
colnames(gridpoint)[1] <- 'LocationName'

gridpoint <- st_transform(
  gridpoint,
  crs = utm_crs
)

# join them
spatial_data <- dplyr::left_join(
  gridpoint,
  tree_imperv,
  by = "LocationName"
) %>% 
  left_join(
    .,
    housing,
    by = "LocationName"
  )

# read in county layers
county <- sf::st_read(
  "D:/GIS/illinois_county",
  layer = "IL_BNDY_County_Py"
)

# Query down to the counties we need. We first filter down to 
#  county name and then select the single COUNTY_NAM (county name)
#  column from the data.
county <- county %>%
  filter(
    .,
    COUNTY_NAM %in% c(
      "MCHENRY",
      "LAKE",
      "KANE",
      "DUPAGE",
      "COOK",
      "KENDALL",
      "WILL"
    )
  ) %>% 
  select(
    COUNTY_NAM
  )

county <- sf::st_transform(
  county,
  crs = utm_crs
)

# remove points outside of our study area
my_raster <- sf::st_intersection(
  spatial_data,
  county
)

# and reduce the county lines to the study area

county_crop <-  sf::st_crop(
  county,
  spatial_data
) %>% 
  sf::st_buffer(., 0)

# create the urbanization metrics
source("./R/prep_data.R")

# read in the results as well
base_results <- read.csv(
  "./results/parameter_estimates.csv",
  stringsAsFactors = FALSE,
  row.names = 1
)

# set up the rotation correctly
if(urb_cov$rotation[1,1] < 0){
urb_cov$rotation[,1] <- urb_cov$rotation[,1] * -1
}

# get the covariates
my_covars <- data.frame(
  house = my_raster$HU10,
  tree = my_raster$tree,
  imp = my_raster$imp
)

# scale the raster data the same way we did in our model
#  The 'urb' df has our unscaled values, and has the 
#  same column names as my_covars
to_scale <- colnames(
  my_covars
)

for(i in 1:3){
  my_covars[,to_scale[i]] <-
    (my_covars[,to_scale[i]] - mean(urb[,to_scale[i]])) /
    sd(urb[,to_scale[i]])
}

# Now we need to apply the PCA to this. This just requires a bit
#  of matrix math (i.e., multiply the rotations by the values). Again
#  everything is in the same order.
urbPCA1 <- as.matrix(my_covars) %*% urb_cov$rotation[,1]
urbPCA2 <- as.matrix(my_covars) %*% urb_cov$rotation[,2]
urbPCA1x2 <- urbPCA1 * urbPCA2

# coyote occupancy. 
coy_occ <- cbind(
  1,
  urbPCA1,
  urbPCA2,
  urbPCA1x2
) %*% base_results[1:4,2] 

# convert to probability
coy_occ <- plogis(
  coy_occ
)

# add this to the raster
my_raster$coyote_occupancy <- coy_occ

# calculate the likelihood of mange now
coy_man <- cbind(
  1,
  urbPCA1,
  urbPCA2,
  urbPCA1x2
) %*% base_results[7:10,2] 

coy_man <- plogis(coy_man)

my_raster$coyote_w_mange <- coy_occ * coy_man

# Look at observed coyote occupancy
site_coords <- read.csv(
  "./data/raw_site_coords.csv",
  stringsAsFactors = FALSE
)

# coyote present, plus an indicator if mange was detected.
coy_p <- coy %>%
  group_by(
    site,
    surveyid
) %>% 
  summarise(
    mange = any(Mange_signs_present > 0)
)

# sdet is the number of seasons detected
# mdet is the number of seasons that mange was detected
coy_p <- coy_p %>% 
  group_by(
    site
  ) %>% 
  summarise(
    sdet = length(surveyid),
    mdet = sum(mange)
  )

# join coy_p with the urbanization covariates
coy_occstat <- left_join(
  site_coords,
  coy_p,
  by = "site"
)

# remove sites not sampled
coy_occstat <- coy_occstat[which(coy_occstat$site %in% urb$site),]

# The color and pch for each point
coy_occstat$col <- "gray20"
coy_occstat$col[is.na(coy_occstat$sdet)] <- "white"

coy_occstat$cex <- 1
coy_occstat$cex[between(coy_occstat$sdet, 6, 10)] <- 1.5
coy_occstat$cex[between(coy_occstat$sdet, 10, 15)] <- 2
coy_occstat <- st_as_sf(
  coy_occstat,
  coords = 2:3,
  crs = utm_crs
)

# do the same for just mangy coyote
coy_mange <- inner_join(
  coy_p,
  site_coords,
  by = 'site'
)

# make new columns that denote shape and color of these points
coy_mange$col <- ifelse(
  coy_mange$mdet > 0,
  "gray20",
  "white"
)
coy_mange$cex <- 1

coy_mange <- st_as_sf(
  coy_mange,
  coords = 4:5,
  crs = utm_crs
)

# bring in stream data
streams <- read_sf(
  "D:/GIS/IL_streams/IL_Streams_From_100K_DLG_Ln.shp"  
  #"Z:/TransectStudyGIS/Study_design/OtherShapefiles/Streams/IL_Streams_nad83.shp"
)

streams <- st_transform(
  streams,
  crs = utm_crs
)

steams <- st_intersection(
  streams,
  county_crop
)

tiff("coyote_occupancy.tiff", width = 5, height = 5,
     units = "in", res = 600, compression = "lzw")
plot(
  my_raster["coyote_occupancy"],
  pch = 15,
  reset = FALSE,
  cex = 0.5,
  main = NULL,
  key.pos = NULL,
  pal = sf.colors(10, alpha = 0.8)

)

# add the stream layer
plot(
  st_geometry(streams),
  add = TRUE,
  col = scales::alpha("white", 0.5),
  lwd = 1.4
)

plot(
  st_geometry(county_crop),
  col = "NA",
  add = TRUE,
  lwd = 2
)
par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            padin =c(0.65, 0.01), lwd = 2, label.cex = 1)
addnortharrow(pos = "topright", padin = c(1.1,0.19), scale = 0.7)

# add points
points(
  st_coordinates(coy_occstat),
  pch = 21,
  bg = scales::alpha(coy_occstat$col, 0.7),
  cex = coy_occstat$cex
)
dev.off()

# just to get the scale bar (putting it all together
#  in inkscape).
svg("./plots/coyote_occupancy_with_scale.svg", width = 5, height = 5)
plot(
  my_raster["coyote_occupancy"],
  pch = 15,
  reset = FALSE,
  cex = 0.5,
  main = NULL,
  bty = "n"
)
dev.off()

tiff("./plots/coyote_mange.tiff", width = 5, height = 5,
     units = "in", res = 600, compression = "lzw")

# coyote mange

plot(
  my_raster["coyote_w_mange"],
  pch = 15,
  reset = FALSE,
  cex = 0.5,
  main = NULL,
  key.pos = NULL,
  breaks = seq(0,1,0.1)
)
# add the stream layer
plot(
  st_geometry(streams),
  add = TRUE,
  col = scales::alpha("white", 0.5),
  lwd = 1.4
)

plot(
  st_geometry(county_crop),
  col = "NA",
  add = TRUE,
  lwd = 2
)

par(xpd = NA)
addscalebar(plotepsg = utm_crs, style = "ticks",
            padin =c(0.65, 0.01), lwd = 2, label.cex = 1)
addnortharrow(pos = "topright", padin = c(1.1,0.19), scale = 0.7)

points(
  st_coordinates(coy_mange),
  pch = 21,
  bg = scales::alpha(coy_mange$col, 0.9),
  cex = coy_mange$cex
)
dev.off()

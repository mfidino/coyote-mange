library(runjags)
library(coda)

source("./R/prep_data.R")

mout <- readRDS("./results/coyote_mcmc_autologistic.RDS")
mout <- as.mcmc.list(mout)

ans <- summary(mout)

base_results <- ans$quantiles[,c(1,3,5)]
base_results <- round(base_results, 2)

# move auto-regressive terms up
base_results <- rbind(
  base_results[1:4,],
  base_results[66,],
  base_results[5:10,],
  base_results[67,],
  base_results[11:65,],
  base_results[68:nrow(base_results),]
)

pinfo <- c("psi - intercept",
           "psi - urb 1",
           "psi - urb 2",
           "psi - urb 1x2",
           "psi - AR(1)",
           "rho - intercept",
           "rho - temperature",
           "omega - intercept",
           "omega - urb 1",
           "omega - urb 2",
           "omega - urb 1x2",
           "omega - AR(1)",
           "gamma - intercept",
           "gamma - blur",
           "gamma - in color",
           "gamma - proportion body visible"
           )

row.names(base_results)[1:16] <- pinfo


# determine if the sign of the median and 95% quantiles are
#  the same (i.e., do not overlap zero).
sig_at_95 <- which(
  rowSums(
    sign(
      base_results[1:16,]
    )
  ) %in% c(3, -3)
)

base_results <- cbind(
  base_results,
  0
)

base_results[sig_at_95,4] <- 1
colnames(base_results)[4] <- "95BCI_excludes_zero"    

write.csv(
  base_results,
  "./results/parameter_estimates.csv", 
  quote = FALSE
)

#####################

# make some plots

# psi stuff
tmp1 <- seq(-3.5, 2.5, length.out = 100)
tmp2 <- seq(-1.5, 2, length.out = 100)

X <- expand.grid(
  tmp1,
  tmp2
)

# make the design matrix
X <- as.matrix(
  cbind(
    1,
    X,
    X[,1] * X[,2]
    )
)

# get just the psi and omega results
psi <- base_results[1:4,2]
psi_theta <- sum(base_results[c(1,5), 2])
ome <- base_results[8:11, 2]
ome_theta <- sum(base_results[c(8,12),2])
pp <- plogis(X %*% psi) / (plogis(X %*% psi) + (1 - plogis(X %*% psi + psi_theta)))
oo <- plogis(X %*% ome) / (plogis(X %*% ome) + (1 - plogis(X %*% ome + ome_theta)))



df <- cbind(X, pp)
dfo <- cbind(X, oo)

# a matrix of the mange results
res_mat <- matrix(
  dfo[,5] * df[,5],
  nrow = length(tmp1),
  ncol =length(tmp2)
)

# determine where...
#  1) coyotes were detected and 
#  2) coyotes with mange were detected

# coyote present, plus an indicator if mange was detected.
coy_p <- coy %>% 
  dplyr::group_by(site, surveyid) %>% 
  dplyr::summarise(
    mange = any(Mange_signs_present > 0)
  )

# sdet is the number of seasons detected
# mdet is the number of seasons that mange was detected
coy_p <- coy_p %>% 
  dplyr::group_by(site) %>% 
  dplyr::summarise(
    sdet = length(surveyid),
    mdet = sum(mange)
  )

y_det[which(y_det$station == 'S03-PSP1'),]

# join coy_p with the urbanization covariates
coy_p <- dplyr::inner_join(
  coy_p,
  y_det[,c('station', 'urb1', 'urb2')],
  by = c('site' = 'station')
) %>% 
  dplyr::distinct(.)



# make new columns that denote shape and color of these points
coy_p$color <- ifelse(
  coy_p$mdet > 0,
  "gray20",
  "gray"
)
coy_p$pch <- ifelse(
  coy_p$mdet > 0,
  24,
  21
)

purple_pal <- colorRampPalette(
  c('white', '#C083F8', '#8F32CE')
)

# plot for coyote with mange
tiff(
  './plots/coyote_mange_plot.tiff',
  height = 5,
  width = 6, 
  units = 'in',
  res = 600,
  compression = "lzw"
)
# Bracket put here so it runs through all the plotting
#  code for this figure at once (in case you are running
#  through this line by line).
{par(mar= c(4,5.5,1,1))
plot(
  1~1, 
  type = 'n', 
  ylim = range(tmp2),
  xlim = range(tmp1), 
  bty = 'n',
  xaxt = 'n', 
  yaxt = 'n', 
  xlab = "", 
  ylab = "", 
  xaxs = 'i', 
  yaxs = 'i'
)
image(
  tmp1,
  tmp2,
  res_mat,
  col = purple_pal(100),
  xlab = "URB1",
  ylab = "URB2",
  las = 1,
  xaxt = 'n',
  yaxt = 'n',
  add = TRUE
)
axis(
  1,
  seq(min(tmp1), max(tmp1), 1),
  labels = FALSE,
  tck = -0.035
)
axis(
  1,
  seq(min(tmp1), max(tmp1), 0.5),
  labels = FALSE,
  tck = -0.035/2
)
mtext(
  text = sprintf(
    "%.1f",
    seq(min(tmp1), max(tmp1), 1)
  ),
  1, 
  line = 0.9,
  at = seq(min(tmp1), max(tmp1), 1),
  cex = 1.2
)
axis(
  2,
  seq(min(tmp2), max(tmp2), 0.5),
  labels = FALSE,
  tck = -0.035,
  cex = 1.1
)
axis(
  2,
  seq(min(tmp2), max(tmp2), 0.25),
  labels = FALSE,
  tck = -0.035/2,
  cex = 1.1
)
mtext(
  text = sprintf(
    "%.1f",
    seq(min(tmp2), max(tmp2), 0.5)
  ),
  2,
  line = 0.9,
  at = seq(min(tmp2), max(tmp2), 0.5),
  cex = 1.2,
  las = 1
)
mtext(
  text = "URB1",
  1,
  line = 2.5,
  at = mean(tmp1),
  cex = 1.5
)
mtext(
  text = "URB2",
  2, 
  line = 3.3,
  at = mean(tmp2),
  cex = 1.5
)
contour(
  tmp1,
  tmp2,
  res_mat,
  add = TRUE,
  nlevels = 5,
  lwd = 2,
  labcex = 1.2
)
points(
  coy_p$urb2 ~ coy_p$urb1,
  bg = coy_p$color,
  pch = coy_p$pch,
  cex = 1.2
)
legend(
  'bottomleft',
  bg = 'white',
  pt.bg = c('gray20', 'gray'),
  pch = c(24,21), 
  legend = c('with mange', 'without mange'),
  title = 'Coyote detected...',
  pt.cex = 1.2, cex = 1.1, y.intersp = 1.5
)
}
dev.off()

# coyotes in general plot

# a matrix of the coyote results
res_mat <- matrix(
  df[,5],
  nrow = length(tmp1),
  ncol =length(tmp2)
)

coy_occ <- y_det[complete.cases(y_det),] %>% 
  dplyr::group_by(station) %>% 
  dplyr::summarise(coy = as.numeric(any(y >0)),
                   urb1 = unique(urb1),
                   urb2 = unique(urb2))
                  

# make new columns that denote shape and color of these points
coy_occ$color <- ifelse(
  coy_occ$coy > 0,
  "gray20",
  "gray"
)
coy_occ$pch <- ifelse(
  coy_occ$coy > 0,
  23,
  22
)

# coyote, mangy or otherwise, plot
tiff(
  './plots/coyote_occupancy_plot.tiff',
  height = 5,
  width = 6, 
  units = 'in',
  res = 600,
  compression = "lzw"
)
{par(
  mar= c(4,5.5,1,1)
)
plot(
  1~1, 
  type = 'n',
  ylim = range(tmp2),
  xlim = range(tmp1),
  bty = 'n',
  xaxt = 'n',
  yaxt = 'n',
  xlab = "",
  ylab = "",
  xaxs = 'i',
  yaxs = 'i'
)
image(
  tmp1,
  tmp2,
  res_mat,
  col = purple_pal(100),
  xlab = "URB1",
  ylab = "URB2",
  las = 1,
  xaxt = 'n',
  yaxt = 'n',
  add = TRUE
)
axis(
  1,
  seq(min(tmp1), max(tmp1), 1),
  labels = FALSE,
  tck = -0.035
)
axis(
  1,
  seq(min(tmp1), max(tmp1), 0.5),
  labels = FALSE,
  tck = -0.035/2
)
mtext(
  text = sprintf(
    "%.1f",
    seq(min(tmp1), max(tmp1), 1)
  ),
  1,
  line = 0.9,
  at = seq(min(tmp1), max(tmp1), 1),
  cex = 1.2
)
axis(
  2,
  seq(min(tmp2), max(tmp2), 0.5),
  labels = FALSE,
  tck = -0.035,
  cex = 1.1
)
axis(
  2,
  seq(min(tmp2), max(tmp2), 0.25),
  labels = FALSE,
  tck = -0.035/2,
  cex = 1.1
)
mtext(
  text = sprintf(
    "%.1f",
    seq(min(tmp2), max(tmp2), 0.5)
  ),
  2,
  line = 0.9,
  at = seq(min(tmp2), max(tmp2), 0.5),
  cex = 1.2,
  las = 1
)
mtext(
  text = "URB1",
  1,
  line = 2.5,
  at = mean(tmp1),
  cex = 1.5
)
mtext(
  text = "URB2",
  2,
  line = 3.3,
  at = mean(tmp2),
  cex = 1.5
)
contour(
  tmp1,
  tmp2,
  res_mat,
  add = TRUE,
  nlevels = 5,
  lwd = 2, 
  labcex = 1.2
)
points(
  coy_occ$urb2 ~ coy_occ$urb1,
  bg = coy_occ$color,
  pch = coy_occ$pch, 
  cex = 1.2
)
legend(
  'bottomleft',
  bg = 'white',
  pt.bg = c('gray20', 'gray'),
  pch = c(23,22), 
  legend = c(
    'Coyote\ndetected',
    'Coyote not\ndetected'
  ),
  pt.cex = 1.2,
  cex = 1.1,
  y.intersp=1.5
)
}
dev.off()

# Coyote with mange over time
n_mange <- base_results[grep("mange", row.names(base_results)),]


{tiff('./plots/site_mange.tiff', height = 5, width = 7, res = 600, units = 'in',
     compression = "lzw")
par(mar= c(4,5.5,1,1))
plot(n_mange[,2], type = 'p', ylim = range(n_mange), bty = 'l', pch = 16,
     ylab = "", xlab = "", xaxt = 'n', yaxt = 'n', cex = 1.2)
for(i in 1:16){
  lines(x = rep(i, 2), y = n_mange[i,c(1,3)])
}
lines(n_mange[,2] ~c(1:16) , lty = 2)
se <- c("SP", "SU", "FA", "WI")
axis(1, at= c(1:16), labels = F, tck = -0.035/2)
mtext("Year:", 1, line = 0.4, at = -0.9, cex = 1.2)
mtext("Season:", 1, line = 1.5, at = -1.28, cex = 1.2)
mtext(text = c(2010:2014), 1, line = 0.4, at = c(1,4,8,12, 16), cex = 1)
mtext(text = rep(se, 4),1, line =1.5, at = c(1:16), cex = 1)
axis(2, at = seq(0,40, 10), labels = F, tck = -0.035/2)
axis(2, at = seq(0,40, 5), labels = F, tck = -0.035/4)
mtext(text = seq(0,40,10),2, las = 1, at = seq(0,40,10), line = 0.8)

coy_sea <- coy %>% group_by(site, season) %>% 
  summarise(mange = any(Mange_signs_present>0)) %>% 
  ungroup(.) %>% group_by(season) %>% 
  summarise(mange = sum(mange))

points(coy_sea, pch = 21, bg = "gray", cex = 1.2)
mtext(text = "Sites with mange", 2, at = 20, line = 3, cex = 1.4)
legend('topright', pt.bg = c('black', 'gray'), pch = c(21), 
       legend = c('Estimated', 'Observed'),
       pt.cex = 1.2, cex = 1.1, bty = "n")
dev.off()
}



# Coyote over time
n_coyote <- base_results[grep("coyote", row.names(base_results)),]
y_det$season <- factor(
  y_det$season,
  levels = levels(coy$season)
)

{tiff('./plots/site_coyote.tiff', height = 5, width = 7, res = 600, units = 'in',
     compression = "lzw")
par(mar= c(4,5.5,1,1))
plot(n_coyote[,2], type = 'p', ylim = c(0,80), bty = 'l', pch = 16,
     ylab = "", xlab = "", xaxt = 'n', yaxt = 'n', cex = 1.2)
for(i in 1:16){
  lines(x = rep(i, 2), y = n_coyote[i,c(1,3)])
}
lines(n_coyote[,2] ~c(1:16) , lty = 2)
se <- c("SP", "SU", "FA", "WI")
axis(1, at= c(1:16), labels = F, tck = -0.035/2)
mtext("Year:", 1, line = 0.4, at = -0.9, cex = 1.2)
mtext("Season:", 1, line = 1.5, at = -1.28, cex = 1.2)
mtext(text = c(2010:2014), 1, line = 0.4, at = c(1,4,8,12, 16), cex = 1)
mtext(text = rep(se, 4),1, line =1.5, at = c(1:16), cex = 1)
axis(2, at = seq(0,80, 10), labels = F, tck = -0.035/2)
axis(2, at = seq(0,80, 5), labels = F, tck = -0.035/4)
mtext(text = seq(0,80,10),2, las = 1, at = seq(0,80,10), line = 0.8)


coy_sea <- y_det %>% group_by(station, season) %>% 
  summarise(coy = any(y>0)) %>% 
  ungroup(.) %>% group_by(season) %>% 
  summarise(coy = sum(coy, na.rm = TRUE))

points(coy_sea, pch = 21, bg = "gray", cex = 1.2)
mtext(text = "Sites with coyote", 2, at = 40, line = 3, cex = 1.4)
legend('bottomright', pt.bg = c('black', 'gray'), pch = c(21), 
       legend = c('Estimated', 'Observed'),
       pt.cex = 1.2, cex = 1.1, bty = "n")
dev.off()
}


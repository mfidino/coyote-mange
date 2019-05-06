
library(dplyr)

# Coyote mange detection in camera trap photos

data <- read.csv("coydata.csv")

# little bit of cleaning
data$Season <- trimws(data$Season)
data$Site <- trimws(data$Site)


# just grabbing the stuff I need to model, keeping the OG data in the 
#  same format just in case I need to query anything.

#  I also logged the blur data because the spread in it sucks.
#   On top of this I discretized it, the average is about
#   560 so I just went down to 500. So these are like
#   mostly clear images?
to_mod <- data.frame(mange = data$Mange_signs_present,
                     color = data$In_color,
                     blur = scale(log(data$blur)),
                     blur500 = as.numeric(log(data$blur) > 6.11),
                     propvis = scale(data$propbodyvis))
                     

mod1 <- glm(mange ~ blur + color + propvis,
            family = binomial, data = to_mod)

mod2 <- glm(mange ~ blur500 + color + propvis,
                    family = binomial, data = to_mod)

AIC(mod1, mod2)
# discretized blur is a better model

# interpretation. Color photos are the best thing to get
#  second best is getting most of the coyote in an image
#  third is for a photo that has at least >500 clarity.


true_vals <- seq(0, 0.6, 0.01)

scale_vals <- (true_vals - mean(data$propbodyvis)) / sd(data$propbodyvis)

to_pred <- data.frame(color = rep(c(0,1), each = length(scale_vals)),
                      blur500 = 1,
                      propvis = rep(scale_vals,times = 2))

my_preds <- predict.glm(mod2,newdata = to_pred, type = 'link', se.fit = TRUE )


my_probs <- data.frame(mu = plogis(my_preds$fit),
                    lower = plogis(qnorm(0.025,my_preds$fit, my_preds$se.fit )),
                    upper= plogis(qnorm(0.975,my_preds$fit, my_preds$se.fit )),
                    color = to_pred$color,
                    propvis = rep(true_vals, times = 2))

# assuming we have a clear photo

windows(8,4)
par(mfrow = c(1,2))
# without color plot
tmp <- my_probs[my_probs$color == 0,]

plot(tmp$mu ~ tmp$propvis, type = 'l', ylim = c(0, 1), lwd = 3, bty ='l',
     xlab = "proportion visible b&W",
     ylab = "probability of classifying as mange")
lines(tmp$lower ~ tmp$propvis, lwd = 2, lty = 2)
lines(tmp$upper ~ tmp$propvis, lwd = 2, lty = 2)
rm(tmp)

# with color plot
tmp <- my_probs[my_probs$color == 1,]

plot(tmp$mu ~ tmp$propvis, type = 'l', ylim = c(0, 1), lwd = 3, bty ='l',
     xlab = "proportion visible color",
     ylab = "probability of classifying as mange")
lines(tmp$lower ~ tmp$propvis, lwd = 2, lty = 2)
lines(tmp$upper ~ tmp$propvis, lwd = 2, lty = 2)
rm(tmp)


# predictions based on the photos at a site (note that when we do this
#  in a Bayesian analysis we will propagate that uncertainty, we are
#  currently working with point estimates).

# we can work directly with the to_mod object to do this, which 
#  means we do not need to supply newdata to the predict function

mod_preds <- predict.glm(mod2, type = "response")

# link this to the site info
mod_preds <- data.frame(prob = mod_preds,
                        site = data$Site,
                        season = data$Season,
                        year = data$Year)

# prob detecting mange at least once is
# 1 - product of prob not detecting mange for all images at a site

my_summary <- mod_preds %>% group_by(site, season, year) %>% 
  summarise(prob_mange = 1 - prod(1 - prob))

# maybe a little interesting!
boxplot(my_summary$prob_mange ~ my_summary$season)



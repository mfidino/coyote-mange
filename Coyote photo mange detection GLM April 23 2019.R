
# Coyote mange detection in camera trap photos

data <- read.csv("coydata.csv")

# just grabbing the stuff I need to model, keeping the OG data in the 
#  same format just in case I need to query anything.

#  I also logged the blur data because the spread in it sucks.
#   On top of this I discretized it, the average is about
#   560 so I just went down to 500. So these are like
#   mostly clear images?
to_mod <- data.frame(mange = data$Mange_signs_present,
                     color = data$In_color,
                     blur = scale(log(data$blur)),
                     disc500 = data$blur > 500,
                     propvis = scale(data$propbodyvis))
                     

mod1 <- glm(mange ~ blur + color + propvis,
            family = binomial, data = to_mod)

mod2 <- mod1 <- glm(mange ~ disc500 + color + propvis,
                    family = binomial, data = to_mod)

AIC(mod1, mod2)
# discretized blur is a better model


true_vals <- seq(0, 0.6, 0.01)

scale_vals <- (true_vals - mean(data$propbodyvis)) / sd(data$propbodyvis)

to_pred <- data.frame(color = rep(c(0,1), each = length(scale_vals)),
                      blur = 0,
                      propvis = rep(scale_vals,times = 2))

my_preds <- predict.glm(mod1,newdata = to_pred, type = 'link', se.fit = TRUE )


my_probs <- data.frame(mu = plogis(my_preds$fit),
                    lower = plogis(qnorm(0.025,my_preds$fit, my_preds$se.fit )),
                    upper= plogis(qnorm(0.975,my_preds$fit, my_preds$se.fit )),
                    color = to_pred$color,
                    propvis = rep(true_vals, times = 2))


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



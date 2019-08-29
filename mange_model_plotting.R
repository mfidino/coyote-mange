
library(coda)
library(runjags)
mout <- readRDS("./results/coyote_mcmc_inxs.RDS")
ans <- summary(mout)

mm <- as.matrix(as.mcmc.list(mout), chains = TRUE)



# calculate coyote without mange
newx1 <- seq(-3,3, 0.1)
newx2 <- seq(-3, 3, 0.1)


my_prob <- seq(0,1,0.01)
yo <- cbind(1, expand.grid(newx1, newx2))

psi <- mm[,c(2:4)]
psi <- apply(psi, 2, median)
omega <- mm[,c(7:9)]
omega <- apply(omega, 2, median)
#omega[2] <- 0
psi_pred <- t(plogis(psi %*% t(yo)))
omega_pred <- t(plogis(omega %*% t(yo)))


c_nomange <- psi_pred * (1 - omega_pred)
c_mange <- psi_pred * omega_pred


colfunc <- colorRampPalette(c("blue","white", "red"))
my_col <- data.frame(value = as.character(my_prob), 
                     color = colfunc(length(my_prob)),
                     stringsAsFactors = FALSE)

my_x <- data.frame(value = as.character(round(c_nomange,2)),
                   stringsAsFactors = FALSE)

my_y <- data.frame(value = as.character(round(c_mange, 2)),
                   stringsAsFactors = FALSE)

testx <- left_join(my_x, my_col, by = "value")
testy <- left_join(my_y, my_col, by = "value")


coy_p <- coy %>% group_by(site) %>% 
  summarise( mange = any(Mange_signs_present > 0))

plot(yo[,3] ~ yo[,2], col = testy$color, pch = 15, bty = "l",
     cex = 1.5, xlab = "URB_1_PCA", ylab = "URB_2_PCA")

ah <- 

ah2 <- left_join(urb_mat, coy_p, by = "site")
ah2$col <- "gray"
ah2$col[ah2$mange] <- "black"

points(y = ah2$urb2,x = ah2$urb1, pch = 21, bg = ah2$col)


c_nomange <- t(apply(c_nomange, 2, quantile, probs = c(0.025,0.5,0.975)))
c_mange <- t(apply(c_mange, 2, quantile, probs = c(0.025,0.5,0.975)))

plot(c_nomange[,2] ~ newx[,2], xlab = "Urbanization",
     ylab = "probability of", ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2)
lines(c_nomange[,1] ~ newx[,2], lty = 2)
lines(c_nomange[,3] ~ newx[,2], lty = 2)

lines(c_mange[,2] ~ newx[,2], col = "red", lwd = 2)
lines(c_mange[,1] ~ newx[,2], col = "red", lty = 2)
lines(c_mange[,3] ~ newx[,2], col = "red", lty = 2)

legend('topright', legend = c('Coyote without mange',
                              'Coyote with mange'),
       col = c('black', 'red'), lwd = 2)

# omega pred
omega_2 <- t(apply(omega_pred, 2, quantile, probs = c(0.025,0.5,0.975)))

plot(omega_2[,2] ~ newx[,2], xlab = "Urbanization",
     ylab = "probability of a coyote having mange", 
     ylim = c(0,1), type = 'l', bty = 'l',
     lwd = 2)
lines(omega_2[,1] ~ newx[,2], lty = 2)
lines(omega_2[,3] ~ newx[,2], lty = 2)

esties <- ans[1:28,1:3]

pr_coy <- ans[1:6,2]

omega_rans <- mm[,c(43:58)] 
for(i in 1:nrow(omega_rans)){
  omega_rans[i,] <- mm[i,7] + omega_rans[i,]
}

omega_prob <- apply(plogis(omega_rans), 2, quantile,
                    probs = c(0.025,0.5,0.975))
plot(omega_prob[2,] ~ c(1:16), bty = 'l', ylim = c(0,0.75),
     pch = 20, cex = 2, ylab= 'probability of mange')

for(i in 1:16){
  lines(x = rep(i, 2),
        y = c(omega_prob[1,i], omega_prob[3,i]))
}


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

n_mange <- mm[,grep('n_mange', colnames(mm))]

n_mange <- apply(n_mange, 2, quantile, probs = c(0.025,0.5,0.975))


plot(n_mange[2,], type = 'p', ylim = range(n_mange), bty = 'l', pch = 16,
     ylab = "sites with mange", xlab = "season")
for(i in 1:16){
  lines(x = rep(i, 2), y = n_mange[-2,i])
}

coy_sea <- coy %>% group_by(site, season) %>% 
  summarise(mange = any(Mange_signs_present>0)) %>% 
  ungroup(.) %>% group_by(season) %>% 
  summarise(mange = sum(mange))

points(coy_sea, pch = 16, col = "red")
library(runjags)
library(coda)
mout <- readRDS("./results/final_mange_model.rds")
mout <- as.mcmc.list(mout)

ans <- summary(mout)

base_results <- ans$quantiles[,c(1,3,5)]

pinfo <- c("psi - intercept",
           "psi - urb 1",
           "psi - urb 2",
           "rho - intercept",
           "rho - temperature",
           "omega - intercept",
           "omega - urb 1",
           "omega - urb 2",
           "gamma - intercept",
           "gamma - blur",
           "gamma - in color",
           "gamma - proportion body visible")
row.names(base_results)[1:12] <- pinfo

hm <- which(rowSums(sign(base_results[1:12,])) %in% c(3, -3))

base_results <- cbind(base_results, 0)
base_results[hm,4] <- 1
colnames(base_results)[4] <- "Significant"    

write.csv(base_results, "./results/parameter_estimates.csv", 
          quote = FALSE)

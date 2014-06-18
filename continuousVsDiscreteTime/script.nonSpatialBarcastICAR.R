##
## This model seems to be working well, maybe reduce the degrees of freedom in the Wishart matrix and run for longer time to get convergence
##

rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

library(dplR)
library(mvtnorm)
library(MASS)
suppressMessages(library(MCMCpack))
source('~/barcast/functions/rMVN.R')
source('~/barcast/continuousVsDiscreteTime//mcmc.nonSpatialBarcastICAR.R')

##
## load data
##

pdsi <- read.table('~/barcast/datafiles/Hudson/namerica-pdsi-act.txt', sep = "", header = TRUE)
Y.crn <- read.crn('~/barcast/datafiles/Hudson/NYCpluvial2010a.allCRNsFinalRmHeaders', header = TRUE)
WP <- Y.crn[, 1:32] - 1
WP[is.na(WP)] <- 0
WP <- as.matrix(WP)

t <- dim(WP)[1]
t.o <- 452:555 # years PDSI is observed

WI <- rep(0, length = t)
WI[t.o] <- pdsi$X261

p <- dim(WP)[2]
W <- cbind(WI, WP)
n <- dim(W)[2]
HI <- rep(0, t)
HI[WI != 0] <- 1
HP <- matrix(0, nrow = t, ncol = p)
HP[WP != 0] <- 1

##
## setup priors
##

n.mcmc <- 50000 #50000
n.burn <- floor(n.mcmc / 5)

phi.squared.prior <- 400
delta.0 <- rep(1, p)
delta.1 <- rep(1, p)
num.neighbors <- 1

params <- list(n.mcmc = n.mcmc, phi.squared.prior = phi.squared.prior, delta.0 = delta.0, delta.1 = delta.1, num.neighbors = num.neighbors)

##
## fit MCMC
##

start <- Sys.time()
# Rprof(filename = '~/barcast/nonSpatialBarcastWishartTakeTwo/Rprof.out', line.profiling = TRUE)
out <- mcmc(WI, WP, HI, HP, params, method = 'continuous')
# Rprof(NULL)
# summaryRprof('~/barcast/nonSpatialBarcastWishartTakeTwo/Rprof.out', lines = "show")
finish <- Sys.time() - start
finish 

# jpeg(file = '~/barcast/plots/DiscreteTimeJune18.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')
layout(matrix(1:9, 3, 2))
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('blue', alpha.f = 0.5))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025), col = adjustcolor('red', alpha.f = 0.25))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975), col = adjustcolor('red', alpha.f = 0.25))

matplot(WI[t.o], type = 'l', col = adjustcolor('blue', alpha.f = 0.5))
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean)[t.o], type = 'l', add = TRUE)
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025)[t.o], col = adjustcolor('red', alpha.f = 0.25))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975)[t.o], col = adjustcolor('red', alpha.f = 0.25))

matplot(out$sigma.squared.save[n.burn:n.mcmc], type = 'l', main = paste('sigma.squared = ', round(mean(out$sigma.squared.save[n.burn:n.mcmc]), digits = 4)))
matplot(out$tau.squared.I.save[n.burn:n.mcmc], type = 'l', main = paste('tau.squared.I = ', round(mean(out$tau.squared.I.save[n.burn:n.mcmc]), digits = 4)))
matplot(out$tau.squared.P.save[n.burn:n.mcmc], type = 'l', main = paste('tau.squared.P = ', round(mean(out$tau.squared.P.save[n.burn:n.mcmc]), digits = 4)))

matplot(apply(out$trend.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('red', alpha.f = 0.5), ylim = c(-0.5, 0.5))
lines(apply(out$trend.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025), type = 'l', col = adjustcolor('red', alpha.f = 0.25))
lines(apply(out$trend.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975), type = 'l', col = adjustcolor('red', alpha.f = 0.25))
# dev.off()


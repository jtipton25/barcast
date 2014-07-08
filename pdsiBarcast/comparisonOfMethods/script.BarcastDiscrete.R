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
source('~/barcast/continuousVsDiscreteTime/mcmc.nonSpatialBarcastICAR.R')
source('~/barcast/plots/make.output.plot.CVsD.R')

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

n.mcmc <- 25000
n.burn <- floor(n.mcmc / 5)

## prior for tau$2_P
alpha.P <- 1
beta.P <- 3
curve(dinvgamma(x, alpha.P, beta.P), 0, 10)

## prior for tau^2_I
alpha.I <- 10
beta.I <- 1
curve(dinvgamma(x, alpha.I, beta.I), 0, 10)

## prior for sigma^2
alpha.sigma <- 2
beta.sigma <- 2
curve(dinvgamma(x, alpha.sigma, beta.sigma), 0, 10)

phi.squared.lower <- 30
phi.squared.upper <- 50
N.phi <- 100
phi.tune <- 1 # potential number of steps to take for proposal of discrete uniform
delta.0 <- rep(1, p)
delta.1 <- rep(1, p)
num.neighbors <- 10
num.truncate <- 100

params <- list(n.mcmc = n.mcmc, phi.squared.lower = phi.squared.lower, phi.squared.upper = phi.squared.upper, N.phi = N.phi, phi.tune = phi.tune, delta.0 = delta.0, delta.1 = delta.1, num.neighbors = num.neighbors, num.truncate = num.truncate, alpha.P = alpha.P, beta.P = beta.P, alpha.I = alpha.I, beta.I = beta.I, alpha.sigma = alpha.sigma, beta.sigma = beta.sigma)

##
## fit MCMC
##

out <- mcmc(WI, WP, HI, HP, params, 'discrete')

make.output.plot(out, method = 'discrete')

# jpeg(file = '~/barcast/plots/hudsonValleyPDSIdiscrete.jpeg', width = 9, height = 6, quality = 100, res  = 600, units = 'in')  
layout(matrix(1:2, nrow = 2))
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('blue', alpha.f = 0.5), main = 'Reconstruction of PDSI', ylab = 'PDSI', xlab = 'Year')
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025), col = adjustcolor('red', alpha.f = 0.25))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975), col = adjustcolor('red', alpha.f = 0.25))
abline(h = 0)

matplot(WI[t.o], type = 'l', col = adjustcolor('blue', alpha.f = 0.5), main = 'Calibration Interval', ylab = 'PDSI', xlab = 'Year')
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean)[t.o], type = 'l', add = TRUE)
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025)[t.o], col = adjustcolor('red', alpha.f = 0.25))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975)[t.o], col = adjustcolor('red', alpha.f = 0.25))
abline(h = 0)
# dev.off()


# start <- Sys.time()
# out <- mcmc(WI, WP, HI, HP, params, method = 'discrete')
# out <- mcmc(WI, WP, HI, HP, params, method = 'continuous')
# finish <- Sys.time() - start
# finish 

# make.output.plot(out, method = 'continuous', file = '~/barcast/plots/ContinuousFitJune19.jpeg')
# make.output.plot(out, method = 'continuous')

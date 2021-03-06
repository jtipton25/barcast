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

n.mcmc <- 5000#25000
n.burn <- floor(n.mcmc / 5)

## prior for tau$2_P
alpha.P <- 1
beta.P <- 3
curve(dinvgamma(x, alpha.P, beta.P), 0, 10)

## prior for tau^2_I
alpha.I <- 100
beta.I <- 1
curve(dinvgamma(x, alpha.I, beta.I), 0, 10)

## prior for sigma^2
alpha.sigma <- 2
beta.sigma <- 2
curve(dinvgamma(x, alpha.sigma, beta.sigma), 0, 10)

phi.lower <- 30
phi.upper <- 50
N.phi <- 100
phi.tune <- 1 # potential number of steps to take for proposal of discrete uniform
delta.0 <- 1 #rep(1, p)
delta.1 <- 1 #rep(1, p)
num.neighbors <- 10
num.truncate <- 100

params <- list(n.mcmc = n.mcmc, phi.lower = phi.lower, phi.upper = phi.upper, N.phi = N.phi, phi.tune = phi.tune, delta.0 = delta.0, delta.1 = delta.1, num.neighbors = num.neighbors, num.truncate = num.truncate, alpha.P = alpha.P, beta.P = beta.P, alpha.I = alpha.I, beta.I = beta.I, alpha.sigma = alpha.sigma, beta.sigma = beta.sigma)

##
## fit MCMC
##

out <- mcmc(WI, WP, HI, HP, params, 'continuous')

make.output.plot(out, method = 'continuous')

# jpeg(file = '~/barcast/plots/hudsonValleyPDSIcontinuous.jpeg', width = 9, height = 6, quality = 100, res  = 600, units = 'in')  
layout(matrix(1:2, nrow = 2))
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('black', alpha.f = 0.5), main = 'Reconstruction of PDSI', ylab = 'PDSI', xlab = 'Year')
lines(WI, , col = adjustcolor('blue', alpha.f = 0.5))
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

# save(out, WI, file = '~/barcast/data/coopMeetingContinuousModelFit10_01_2014.RData')
# load(file = '~/barcast/data/coopMeetingContinuousModelFit.RData')
# load(file = '~/barcast/data/continuousFitOutputThu Jul 17 14:11:41 2014.RData')
# load(file = '~/barcast/data/coopMeetingContinuousModelFit10_01_2014.RData')
year <- 2005 - (556:1)
t <- length(year)


pederson <- read.table(file = '~/Google Drive/NYCDroughtPederson.txt', sep = '', header = TRUE)
pederson$YEAR[1:474]


# jpeg(file = '~/barcast/plots/continuousReconstructionPederson.jpeg', width = 18, height = 6, quality = 100, res  = 600, units = 'in')
matplot(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, median), type = 'l', col = adjustcolor('black', alpha.f = 0.5), main = 'Reconstruction of PDSI', ylab = 'PDSI', xlab = 'Year', xaxt = 'n')
matplot(c(rep(NA, 82), pederson$drought[1:474]), add = TRUE, type = 'l', col = adjustcolor('red', alpha.f = 0.5))
matplot(c(rep(NA, length(WI[WI == 0]) - 1), WI[WI!=0], NA), type = 'l', add = TRUE, col = adjustcolor('blue', alpha.f = 0.5))
# lines(WI, , col = adjustcolor('blue', alpha.f = 0.5))
# lines(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = 0.025), col = adjustcolor('red', alpha.f = 0.25))
# lines(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = 0.975), col = adjustcolor('red', alpha.f = 0.25))

axis(1, at = c(2 + (0:11) * 50), labels = year[c(2 + (0:11) * 50)])
idx <- 0:19
for(i in 1:19){
  polygon(c(1:t, t:1), c(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = (0.5 + (idx[i + 1]) / 40)), rev(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = (0.5 + idx[i] / 40)))), col = adjustcolor('grey80', alpha.f = (1 - idx[i] / 20)), border = NA)  
  polygon(c(1:t, t:1), c(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = (0.5 - (idx[i + 1]) / 40)), rev(apply(out$T.save[, (dim(out$T.save)[2] / 5) : dim(out$T.save)[2]], 1, quantile, prob = (0.5 - idx[i] / 40)))), col = adjustcolor('grey80', alpha.f = (1 - idx[i] / 20)), border = NA)  
}
# abline(h = 0, col = adjustcolor('black', alpha.f = 0.5))
# dev.off()



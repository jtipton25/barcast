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
source('~/barcast/continuousVsDiscreteTime/mcmc.continuous.one.beta.R')
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

##
## Calibration and validation
##

t.val <- 500:555
t.cal <- 452:499
WI.t.val <- WI[t.val]
WI[t.val] <- 0

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


## prior for tau$2_P
alpha.P <- 0.1 #3
beta.P <- 0.1 #10
curve(dinvgamma(x, alpha.P, beta.P), 0, 10)

## prior for tau^2_I
alpha.I <- 0.1 #200
beta.I <- 0.1 #2
curve(dinvgamma(x, alpha.I, beta.I), 0, 10)

## prior for sigma^2
alpha.sigma <- 0.1
beta.sigma <- 0.1
curve(dinvgamma(x, alpha.sigma, beta.sigma), 0, 10)


phi.tune <- 1 # potential number of steps to take for proposal of discrete uniform
delta.0 <- 1 # rep(1, p)
delta.1 <- 1 # rep(1, p)
num.neighbors <- 10
num.truncate <- t

phi.lower <- 1
phi.upper <- 1000
N.phi <- 1000
n.mcmc <- 10000
n.burn <- floor(n.mcmc / 5)

params <- list(n.mcmc = n.mcmc, phi.lower = phi.lower, phi.upper = phi.upper, N.phi = N.phi, phi.tune = phi.tune, delta.0 = delta.0, delta.1 = delta.1, num.neighbors = num.neighbors, num.truncate = num.truncate, alpha.P = alpha.P, beta.P = beta.P, alpha.I = alpha.I, beta.I = beta.I, alpha.sigma = alpha.sigma, beta.sigma = beta.sigma)

##
## fit MCMC
##

start <- Sys.time()
out <- mcmc.cont(WI, WP, HI, HP, params)
finish <- Sys.time() - start
finish

make.output.plot(out, method = 'continuous')

# jpeg(file = '~/barcast/plots/ContinuouscomparisonPlot.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')  
layout(matrix(1:3, nrow = 3, ncol = 1))
matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean)[ - 1], type = 'l', main = paste('fitted PDSI, runtime for 10000 iterations = ', round(finish, digits = 2), " minutes", sep = ""), ylab = 'PDSI', ylim = c(-4, 4))
matplot(apply(out$trend.save[, n.burn:n.mcmc], 1, mean), type = 'l', add = TRUE, col = 'dark green')
abline(h = 0, col = 'blue')
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025)[ - 1], col = adjustcolor('red', alpha = 0.25))
lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975)[ - 1], col = adjustcolor('red', alpha = 0.25))
lines(WI, col = adjustcolor('blue', alpha = 0.5))
MSPE.val <- (apply(out$T.save[, n.burn:n.mcmc], 1, mean)[t.val] - WI.t.val)^2
matplot(MSPE.val, type = 'l', main = paste('MSPE error over validation', round(mean(MSPE.val), digits = 4)))
matplot(WI.t.val, type = 'l', col = 'blue', main = 'Plot of fitted vs. held out PDSI')
lines(apply(out$T.save[, n.burn:n.mcmc], 1, mean)[t.val], col = 'black')


# dev.off()

# save(out, file = paste('~/barcast/data/continuousFitOutput', date(), '.RData', sep = ''))

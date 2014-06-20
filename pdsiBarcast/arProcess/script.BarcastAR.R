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
source('~/barcast/plots/make.output.plot.wishart.R')
source('~/barcast/pdsiBarcast/arProcess/mcmc.nonSpatialBarcast.R')

##
## load data
##

pdsi <- read.table('~/barcast/datafiles/Hudson/namerica-pdsi-act.txt', sep = "", header = TRUE)
Y.crn <- read.crn('~/barcast/datafiles/Hudson/NYCpluvial2010a.allCRNsFinalRmHeaders', header = TRUE)
# n.trees <- read.table('~/treeRing/datafiles/Hudson/num.trees.txt')
WP <- Y.crn[, 1:32] - 1
# Y[idx.na] <- Y[idx.na] - 1
years <- 1449:2004 # unobserved years include 1449:1899 and 2004

t <- dim(WP)[1]
# n.tree.mat <- matrix(NA, nrow = t, ncol = n)
# n.tree.mat[idx.na] <- n.trees$V1
t.u <- c(1:451, 556) #years PDSI is unobserved
t.o <- 452:555 # years PDSI is observed
# t.cal <- 500:555 
# t.val <- 452:499
t.val <- 500:555
t.cal <- 452:499
##
## Calibration and validation
##
WI <- rep(NA, length = t)
WI[t.o] <- pdsi$X261
WI.t.val <- WI[t.val]
WI[t.val] <- NA

W <- cbind(WI, WP)
n <- dim(W)[2]
H <- !is.na(W)
H.tmp <- rbind(rep(FALSE, dim(W)[2]), H)

year <- rep(1:t, times = (length(W) / t))

##
## setup priors
##

mu.0.tilde <- 0
sigma.squared.0.tilde <- 1
alpha.alpha <- 0
beta.alpha <- 1
mu.0 <- 0
sigma.squared.0 <- 1
alpha.epsilon <- 1
beta.epsilon  <- 1
alpha.phi <- 0.1
beta.phi <- 0.1
alpha.I <- 0.1
beta.I <- 0.1
alpha.P <- 0.1
beta.P <- 0.1
mu.beta.0 <- 0
sigma.squared.beta.0 <- 8
mu.beta.1 <- 0
sigma.squared.beta.1 <- 8

n.mcmc <- 5000
n.burn <- floor(n.mcmc / 5)

##
## fit MCMC
##

start <- Sys.time()
# Rprof(filename = '~/barcast/nonSpatialBarcastWishartTakeTwo/Rprof.out', line.profiling = TRUE)
out <- mcmc(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1)
# Rprof(NULL)
# summaryRprof('~/barcast/nonSpatialBarcastWishartTakeTwo/Rprof.out', lines = "show")
finish <- Sys.time() - start
finish 

# jpeg(file = '~/barcast/plots/ARcomparisonPlot.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')  
layout(matrix(1:3, nrow = 3, ncol = 1))
matplot(apply(out$X.save[, n.burn:n.mcmc], 1, mean)[ - 1], type = 'l', main = paste('fitted PDSI, runtime for 5000 iterations = ', round(finish, digits = 2), " minutes", sep = ""), ylab = 'PDSI', ylim = c(-4, 4))
abline(h = 0, col = 'blue')
lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025)[ - 1], col = adjustcolor('red', alpha = 0.25))
lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975)[ - 1], col = adjustcolor('red', alpha = 0.25))
lines(WI, col = adjustcolor('blue', alpha = 0.5))
MSPE.val <- (apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.val + 1] - WI.t.val)^2
matplot(MSPE.val, type = 'l', main = paste('MSPE error over validation', round(mean(MSPE.val), digits = 4)))
matplot(WI.t.val, type = 'l', col = 'blue', main = 'Plot of fitted vs. held out PDSI')
lines(apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.val + 1], col = 'black')
# dev.off()

rm(list = ls())

##
## libraries and subroutines
##

library(dplR)
library(mvtnorm)
library(MASS)
source('~/barcast/rMVN.R')
source('~/barcast/plots/make.output.plot.R')
source('~/barcast/nonSpatialBarcast/mcmc.nonSpatialBarcast.R')

##
## load data
##

pdsi <- read.table('~/barcast/datafiles/Hudson/namerica-pdsi-act.txt', sep = "", header = TRUE)
Y.crn <- read.crn('~/barcast/datafiles/Hudson/NYCpluvial2010a.allCRNsFinalRmHeaders', header = TRUE)
# n.trees <- read.table('~/treeRing/datafiles/Hudson/num.trees.txt')
WP <- Y.crn[, 1:32]
# Y[idx.na] <- Y[idx.na] - 1
years <- 1449:2004 # unobserved years include 1449:1899 and 2004

t <- dim(WP)[1]
# n.tree.mat <- matrix(NA, nrow = t, ncol = n)
# n.tree.mat[idx.na] <- n.trees$V1
t.u <- c(1:451, 556) #years PDSI is unobserved
t.o <- 452:555 # years PDSI is observed
WI <- rep(NA, length = t)
WI[t.o] <- pdsi$X261

W <- cbind(WI, WP)
H <- !is.na(W)



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
alpha.sigma.squared <- 0.1
beta.sigma.squared  <- 0.1
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

n.mcmc <- 10000
n.burn <- floor(n.mcmc / 5)

##
## fit MCMC
##

start <- Sys.time()
out <- mcmc(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.sigma.squared, beta.sigma.squared, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1)
finish <- Sys.time() - start
finish 

# x11()
make.output.plot(out)
# make.output.plot(out, model = 'simple', file = '~/treeRing/plots/simple.log.climate.png')

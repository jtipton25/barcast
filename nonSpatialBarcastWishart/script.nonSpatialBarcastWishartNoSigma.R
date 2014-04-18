##
## This model seems to be working well, maybe reduce the degrees of freedom in the Wishart matrix and run for longer time to get convergence
##

rm(list = ls())

##
## libraries and subroutines
##

library(dplR)
library(mvtnorm)
library(MASS)
library(MCMCpack)
source('~/barcast/rMVN.R')
source('~/barcast/plots/make.output.plot.wishart.R')
source('~/barcast/nonSpatialBarcastWishart/mcmc.nonSpatialBarcastWishartNoSigma.R')

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
WI <- rep(NA, length = t)
WI[t.o] <- pdsi$X261

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
alpha.sigma.squared <- 1
beta.sigma.squared  <- 1
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

n^2 - n - t
nu.wish <- n + 2 # Must be between n and n^2 - t

##
## Need to determine how to keep the wishart matrix from growing as the degrees of freedom increase
##

n.mcmc <- 1000
n.burn <- floor(n.mcmc / 5)

##
## fit MCMC
##

start <- Sys.time()
# Rprof(file = '~/barcast/nonSpatialBarcastWishart/Rprof.out', line.profiling = TRUE)
out <- mcmc(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.sigma.squared, beta.sigma.squared, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1, nu.wish)
# Rprof(NULL)
# summaryRprof('~/barcast/nonSpatialBarcastWishart/Rprof.out')
finish <- Sys.time() - start
finish 

# x11()
make.output.plot(out)

# for(i in 1:33){
#   matplot(out$Q.save[i, , n.burn:n.mcmc], type = 'l')
# }

# y <- apply(out$X.save[, 1, n.burn:n.mcmc], 1, mean)
# x <- 1:(t+1)
# abline(lm(y~x), col = 'green')
# summary(lm(y~x))
# make.output.plot(out, file = '~/barcast/plots/wishartCorrelation3.31.2014.jpg')
# make.output.plot(out, model = 'simple', file = '~/treeRing/plots/simple.log.climate.png')
# 
# for(i in 1:32){
#   matplot(apply(out$X.save[, i, n.burn:n.mcmc], 1, mean), type = 'l')
# }
# 
# ##
# ## Plot using the longest tree ring reconstruction
# ##
# 
# layout(1:2)
# matplot(apply(out$X.save[, 13, n.burn:n.mcmc], 1, mean), type = 'l')
# abline(h = 0, col = 'blue')
# lines(apply(out$X.save[, 13, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
# lines(apply(out$X.save[, 13, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
# lines(WI, col = adjustcolor('blue', alpha = 0.5))
# y <- apply(out$X.save[, 13, n.burn:n.mcmc], 1, mean)
# x <- 1:(t+1)
# abline(lm(y~x), col = 'green')
# # 			plot(X, type = 'l', col = 'blue')
# plot(sqrt((WI[t.o] - apply(out$X.save[, 13, n.burn:n.mcmc], 1, mean)[t.o + 1])^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')


X.mat <- matrix(NA, nrow = t + 1, dim(H.tmp)[2])
X.save <- matrix(0, nrow = t + 1, n.mcmc - n.burn)
for(k in (n.burn + 1):n.mcmc){
  if(k %% 100 == 0){
    cat(k, ' ')
  }
  X.mat[H.tmp] <- out$X.save[, , k][H.tmp]
  X.save[, k - n.burn] <- rowMeans(X.mat, na.rm = TRUE)
}

matplot(apply(X.save[-1, ], 1, mean), type = 'l')
abline(h=0, col = 'red')
lines(apply(X.save[-1, ], 1, quantile, prob = 0.025), col = adjustcolor('red', alpha = 0.25))
lines(apply(X.save[-1, ], 1, quantile, prob = 0.975), col = adjustcolor('red', alpha = 0.25))
lines(WI, col = adjustcolor('blue', alpha = 0.5))


##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##

rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

library('dplR')
library('mvtnorm')
library('MASS')
suppressMessages(library('MCMCpack'))
library('statmod')

source('~/barcast/functions/rMVN.R')
source('~/barcast/plots/make.output.plot.wishart.R')
source('~/barcast/nonSpatialBarcastWishartTakeTwo/mcmc.nonSpatialBarcast.R')
## load the ODA mcmc code
source('~/barcast/modelAveraging/mcmc.ODA.no.intercept.R')
## code for plotting the output
source('~/1dSpatialSim/plots/make.output.plot.ci.R')
## code for centering and scaling with variance = SSE / n instead of SSE / (n - 1)
source('~/barcast/functions/scale.predictor.R')

##
## load data
##

pdsi <- read.table('~/barcast/datafiles/Hudson/namerica-pdsi-act.txt', sep = "", header = TRUE)
Y.crn <- read.crn('~/barcast/datafiles/Hudson/NYCpluvial2010a.allCRNsFinalRmHeaders', header = TRUE)
# n.trees <- read.table('~/treeRing/datafiles/Hudson/num.trees.txt')
Y <- as.matrix(Y.crn[, 1:32] - 1)  ## need to incoroporate the number of trees

## unobserved years include 1449:1899 and 2004
years <- 1449:2004 

##years PDSI is unobserved
t.u <- c(1:451, 556) 
## years PDSI is observed
t.o <- 452:555 

## data matrix for years with PDSI observed and for PDSI to be predicteds
Y.o <- Y[t.o, ]
Y.new <- Y[t.u, ]

t <- dim(Y)[1]
p <- dim(Y)[2]
# n.tree.mat <- matrix(NA, nrow = t, ncol = n)
# n.tree.mat[idx.na] <- n.trees$V1

##
## Calibration and validation - not going to consider for now
##

# Y.o <- rep(NA, length = t)
# Y.o[t.o] <- pdsi$X261
Y.o <- pdsi$X261

# year <- rep(1:t, times = (length(W) / t))

H.list <- vector('list', length = t)

for(i in 1:t){
  H.list[[i]] <- which(!is.na(Y[i, ]))
}

##
## setup priors
##

matplot(Y, type = 'l')

D <- vector('list', length = length(t.o))
for(i in 1:t.o){
  D[[i]] <- diag(rep(max(eigen(t(X.o[, t]) %*% X.o)$values), dim(X.o)[2])) + 0.0001
}
X.a <- chol(D - t(X.o) %*% X.o)

X.c <- rbind(X.o, X.a)
t(X.c) %*% X.c


##
## Initialize priors and tuning paramteters
##


alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 50000
lambda <- c(0, rep(1, p))
# lambda <- rep(1, p)
n.burn <- n.mcmc / 5

params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda)

##
## fit mcmc using ODA model
##

out <- mcmc.oda(Y.list = Y.list, X.o = X.o, H.list = H.list, params = params)

## Rao-blackwell estimates

beta.fit <- matrix(nrow = p + 1, ncol = reps / 2)
for(i in 1:(reps / 2)){
  for(j in 1:(p + 1)){
    if(j == 1){
#       beta.fit[j, i] <- apply(out$beta.save[j, i, ], 1, mean)
        beta.fit[j, i] <- mean(out$beta.save[j, i, ])
    } else {
#       beta.fit[j, i] <- apply(
        beta.fit[j, i] <- mean(out$rho.save[j - 1, i, ] * out$delta.save[[i]] / (out$delta.save[[i]] + lambda[2:(p + 1)][j - 1]) * out$beta.save[j - 1, i, ])
        #, 1, mean)  
    }  
  }
}

Y.pred <- matrix(nrow = m, ncol = reps / 2)
for(i in 1:(reps / 2)){
  Y.pred[, i] <- X.o %*% beta.fit[, i]  
#     Y.pred[, i] <- X.pca.scale %*% beta.fit[, i]  
}

matplot(Y.pred, type = 'l')
matplot((Y.pred - X.new)^2, type = 'l')
## mean square prediction error
MSPE <- mean((Y.pred - X.new)^2)
MSPE
# log.score <- mean(out$log.score.save[(n.burn + 1):n.mcmc])


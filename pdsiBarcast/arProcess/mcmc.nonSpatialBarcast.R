##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

## Next: Remove Sigma, Sigma.inv, Sigma.epsilon and Sigma.epsilon.inv
##
##
## Note that in the code X represents the field T in the Barcast Model
mcmc <- function(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.epsilon, beta.epsilon, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1){
  
  ##
  ## libraries and subroutines
  ##
  
  library(truncnorm)
  
  make.Ht <- function(i, beta.1){
    c(HIt[i], (HPt[i, ] * beta.1))
  }
  
  make.Wt <- function(i){
    W[i, ][H.idx[i, ]]
  }
  
  make.WPt <- function(i){
    as.vector(WP[i, ][HP.idx[i, ]])
  }
  
  make.WIt <- function(i){
    if(is.na(WI[i])){
      NULL
    } else {
      as.vector(WI[i][HI.idx[i]])
    }
  }

  ##
  ## initialize variables
  ##
  
  W <- as.matrix(cbind(WI, WP))
  H <- matrix(as.numeric(!is.na(W)), nrow = dim(W)[1], ncol = dim(W)[2])
  H.idx <- !is.na(W)
  HP.idx <- !is.na(WP)
  HI.idx <- !is.na(WI)
  HIt <- H[, 1]
  HPt <- H[, - 1]
  NI <- 1
  NP <- dim(WP)[2]
  NIt <- vector(length = t)
  NPt <- vector(length = t)
  Nt <- vector(length = t)
  for(i in 1:t){
    NIt[i] <- sum(!is.na(WI[i]))
    NPt[i] <- sum(!is.na(WP[i, ]))
    Nt[i] <- NIt[i] + NPt[i]
  }
  MI <- sum(NIt)  
  MP <- sum(NPt)
  W[is.na(W)] <- 0
  W <- as.matrix(W)
  WP[is.na(WP)] <- 0
  WP <- as.matrix(WP)
  WI[is.na(WI)] <- 0
  WI <- as.matrix(WI)
  t <- dim(W)[1]
  n <- dim(W)[2]
  X <- vector(length = t + 1)


  Wt <- lapply(1:t, make.Wt)
  WPt <- lapply(1:t, make.WPt)
  WIt <- lapply(1:t, make.WIt) 
#   HIt <- lapply(1:t, make.HIt)
#   HPt <- lapply(1:t, make.HPt)  
  J <- rep(1, n)
  I <- diag(n)
  Jt <- vector('list', length = t)
  for(i in 1:t){
    Jt[[i]] <- c(0, HPt[i, ])
  }
  
  alpha <- runif(1, alpha.alpha, beta.alpha)
  mu <- rnorm(1, mu.0, sigma.squared.0)
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon)
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  beta.0 <- rnorm(1, mu.beta.0, sigma.squared.beta.0)
  beta.1 <- rnorm(1, mu.beta.1, sigma.squared.beta.1)
  Bt <- vector('list', length = t)
  for(i in 1:t){
    Bt[[i]] <- beta.0 * Jt[[i]]
  }
#   Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
  Sigma.full <- diag(c(rep(tau.squared.I, NI), rep(tau.squared.P, NP)))
  Sigma.full.inv <- diag(c(rep(1 / tau.squared.I, NI), rep(1 / tau.squared.P, NP)))
  Sigma <- diag(c(rep(tau.squared.I, 1), rep(tau.squared.P, n - 1)))
  Sigma.inv <- diag(c(rep(1 / tau.squared.I, 1), rep(1 / tau.squared.P, n - 1)))
  
  ##
  ## set up save variables
  ##
  
  X.save <- matrix(nrow = t + 1, ncol = n.mcmc)
  alpha.save <- vector(length = n.mcmc)
  mu.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc)
  tau.I.save <- vector(length = n.mcmc)
  tau.P.save <- vector(length = n.mcmc)
  beta.0.save <- vector(length = n.mcmc)
  beta.1.save <- vector(length = n.mcmc)
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0) cat(' ', l)
    
    ##
    ## sample X
    ##
    
    ## for t = 0
    A.chol <- sqrt(alpha^2 / sigma.squared.epsilon + 1 / sigma.squared.0.tilde)
    b <- alpha / sigma.squared.epsilon * (X[2] - (1 - alpha) * mu) + mu.0.tilde / sigma.squared.0.tilde
    X[1] <- rMVN(A.chol, b)
    
    # for t = 1,...,T - 1
    for(i in 1:(t - 1)){ # note the offset since t = 0 is the first element of the matrix X
      A.chol <- sqrt(t(c(HIt[i], HPt[i, ] * beta.1)) %*%  Sigma.inv %*% (c(HIt[i], HPt[i, ] * beta.1)) + (alpha^2 + 1) /  sigma.squared.epsilon)
      b <- t(c(HIt[i], HPt[i, ] * beta.1)) %*% Sigma.inv %*% (W[i, ] - Bt[[i]]) + (alpha * X[i] + (1 - alpha) * mu) / sigma.squared.epsilon + alpha * (X[i + 2] - (1 - alpha) * mu) / sigma.squared.epsilon
      X[i + 1] <- rMVN(A.chol, b)
    }
    
    # for t = T
    A.chol <- sqrt(t(H[t, ] * beta.1) %*% Sigma.inv %*% (H[t, ] * beta.1) + 1 / sigma.squared.epsilon)
    b <- t(c(HIt[t], HPt[t, ] * beta.1)) %*% Sigma.inv %*% (W[t, ] - Bt[[t]]) + (alpha * X[t] + (1 - alpha) * mu) / sigma.squared.epsilon
    X[t + 1] <- rMVN(A.chol, b)
    
    ##
    ## sample beta_0
    ##
    
    A.chol <-sqrt(MP / tau.squared.P + 1 / sigma.squared.beta.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- HPt[i, ] %*% (WP[i, ] - beta.1 * HPt[i, ] * X[i + 1])
    }
    b <- sum(tmp) / tau.squared.P + mu.beta.0  / sigma.squared.beta.0
    beta.0 <- rMVN(A.chol, b)
    for(i in 1:t){
      Bt[[i]] <- beta.0 * Jt[[i]]
    }
    
    ##
    ## sample beta_1
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- 
        t(HPt[i, ] * X[i + 1]) %*% (HPt[i, ] * X[i + 1])
    }
    A.chol <-sqrt(sum(tmp) / tau.squared.P + 1 / sigma.squared.beta.1)       
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(HPt[i, ] * X[i + 1]) %*% (WP[i, ] - Bt[[i]][ - 1])
    }
    b <- sum(tmp) / tau.squared.P + mu.beta.1 / sigma.squared.beta.1
    beta.1 <- rMVN(A.chol, b)
#     Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
    
    ##
    ## sample mu
    ##
    
    #     A.chol <- sqrt((1 - alpha)^2 * (MP / tau.squared.P + MI / tau.squared.I) + 1 / sigma.squared.0)
    A.chol <- sqrt(t * (1 - alpha)^2 / sigma.squared.epsilon + 1 / sigma.squared.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <-  (X[i + 1] - alpha * X[i])
    }
    b <- (1 - alpha) * sum(tmp) / sigma.squared.epsilon + mu.0 / sigma.squared.0
    mu <- rMVN(A.chol, b)
    
    ##
    ## sample alpha
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- (X[i] - mu) %*% (X[i] - mu) / sigma.squared.epsilon
    }
    A.inv <- 1 / sum(tmp)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- (X[i] - mu) %*% (X[i + 1] - mu) / sigma.squared.epsilon
    }
    b <- sum(tmp)
    alpha <- rtruncnorm(1, a = 0, b = 1, mean = A.inv * b, sd = sqrt(A.inv))
    
    ##
    ## sample tau.squared.I
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      if(NIt[i] != 0){
        tmp[i] <- (WI[i] - HIt[i] %*% X[i + 1]) %*% (WI[i] - HIt[i] %*% X[i + 1])
      } else {
        tmp[i] <- 0
      }
    }
    tau.squared.I <- 1 / rgamma(1, alpha.I + MI / 2, beta.I + sum(tmp) / 2)
    
    ##
    ## sample tau.squared.P
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(WP[i, ] - beta.1 * HPt[i, ] * X[i + 1] - Bt[[i]][ - 1]) %*% (WP[i, ] - beta.1 * HPt[i, ] * X[i + 1] - Bt[[i]][ - 1])
    }
    tau.squared.P <- 1 / rgamma(1, alpha.P + MP / 2, beta.P + sum(tmp) / 2)
    
    Sigma <- diag(c(rep(tau.squared.I, 1),rep(tau.squared.P, n - 1)))
    Sigma.inv <- diag(c(rep(1 / tau.squared.I, 1),rep(1 / tau.squared.P, n - 1))) 
    
    ##
    ## sample sigma.squared.epsilon
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- (X[i + 1] - alpha * X[i] - (1 - alpha) * mu) %*% (X[i + 1] - alpha * X[i] - (1 - alpha) * mu)
    }
    sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon + t / 2, beta.epsilon + sum(tmp) / 2)
    
    ##
    ## save variables
    ##
    
    X.save[, l] <- X
    alpha.save[l] <- alpha
    mu.save[l] <- mu
    tau.I.save[l] <- tau.squared.I
    tau.P.save[l] <- tau.squared.P
    beta.0.save[l] <- beta.0
    beta.1.save[l] <- beta.1
    sigma.squared.epsilon.save[l] <- sigma.squared.epsilon
  }    
  
  list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save)  
}
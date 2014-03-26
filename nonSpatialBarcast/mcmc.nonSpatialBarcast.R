##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

##
##
##
## Note that in the code X represents the field T in the Barcast Model
mcmc <- function(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.sigma.squared, beta.sigma.squared, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1){
  
  ##
  ## libraries and subroutines
  ##
  
  library(truncnorm)

  make.HIt <- function(i){
    if(NIt[i] == 0){
      return(NULL)
    } else {
      return(c(1, rep(0, n - 1)))
    }
  }

  make.HPt <- function(i){
    tmp <- matrix(0, nrow = NPt[i], ncol = n)
    for(j in 1:NPt[i]){
      tmp[j, 1 + which(HP[i, ] == 1)[j]] <- 1
    }
    return(tmp)
  }
  
  make.Ht <- function(i, beta.1){
    rbind(HIt[[i]], (HPt[[i]] * beta.1))
  }
  
  make.Bt <- function(i, beta.0){
    rbind(HIt[[i]], HPt[[i]]) %*% c(0, rep(beta.0, NP)) 
  }
  
  make.Wt <- function(i){
    W[i, ][H[i,]]
  }
  
  make.WPt <- function(i){
    WP[i, ][HP.idx[i, ]]
  }
  
  make.WIt <- function(i){
    if(is.na(WI[i])){
      NULL
    } else {
      WI[i][HI.idx[i]]
    }
  }
  
  ##
  ## initialize variables
  ##
  
  W <- as.matrix(cbind(WI, WP))
  WP <- as.matrix(WP)
  t <- dim(W)[1]
  n <- dim(W)[2]
#   t.u <- c(1:451, 556) #years PDSI is unobserved
#   t.o <- 452:555 # years PDSI is observed
  X <- matrix(0, nrow = t + 1, ncol = n)
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

  HI.idx <- !is.na(WI)
  HI <- ifelse(HI.idx, 1, 0)
  HP.idx <- !is.na(WP)
  HP <- ifelse(HP.idx, 1, 0)
  H <- !is.na(W)

  Wt <- lapply(1:t, make.Wt)
  WPt <- lapply(1:t, make.WPt)
  WIt <- lapply(1:t, make.WIt)

#   HIt <- vector('list', length = t)
#   for(i in 1:t){
#     if(NIt[i] == 0){
#       HIt[[i]] <- c()
#     } else {
#       HIt[[i]] <- c(1, rep(0, n - 1))
#     }
#   }

  HIt <- lapply(1:t, make.HIt)

#   HPt <- vector('list', length = t)
#   for(i in 1:t){
#     HPt[[i]] <- matrix(0, nrow = NPt[i], ncol = n)
#     for(j in 1:NPt[i]){
#       HPt[[i]][which(HP[i,] == 1)[j]] <- 1
#     }
#   }

  HPt <- lapply(1:t, make.HPt)


#   Ht <- vector('list', length = t)
#   for(i in 1:t){
#     Ht[[i]] <- matrix(0, nrow = Nt[i], ncol = n)
#     for(j in 1:Nt[i]){
#       Ht[[i]][j, which(H[i, ] == 1)[j]] <- 1
#     }
#   }
# 
#   for(i in 1:t){
#     rbind(HIt[[i]], HPt[[i]])
#   }
  
  

  J <- rep(1, n)
  I <- diag(n)
  It <- vector('list', length = t)
  IIt <- vector('list', length = t) 
  IPt <- vector('list', length = t)
  for(i in 1:t){
    IIt[[i]] <- diag(1)
    IPt[[i]] <- diag(NPt[i])
    It[[i]] <- diag(NIt[i] + NPt[i])
  }
  MI <- sum(NIt)  
  MP <- sum(NPt)

  alpha <- runif(1, alpha.alpha, beta.alpha)
  mu <- rnorm(1, mu.0, sigma.squared.0)
  sigma.squared <- 1 / rgamma(1, alpha.sigma.squared, beta.sigma.squared)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  beta.0 <- rnorm(1, mu.beta.0, sigma.squared.beta.0)
  beta.1 <- rnorm(1, mu.beta.1, sigma.squared.beta.1)
  Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
  
  
  Bt <- lapply(1:t, make.Bt, beta.0 = beta.0)

  Sigma.epsilon <- sigma.squared * I
  Sigma.epsilon.inv <- 1 / sigma.squared * I
  Sigma <- vector('list', length = t)
  Sigma.inv <- vector('list', length = t)
  for(i in 1:t){
    Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]),rep(tau.squared.P, NPt[i])))
    Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]),rep(1 / tau.squared.P, NPt[i])))
  }  
#   Sigma <- c(rep(tau.squared.I, NI), rep(tau.squared.P, NP)) * I
#   Sigma.inv <- c(rep(1 / tau.squared.I, NI), rep(1 / tau.squared.P, NP)) * I   


  ##
  ## set up save variables
  ##
  
  X.save <- array(dim = c(t + 1, n, n.mcmc))
  alpha.save <- vector(length = n.mcmc)
  mu.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc)
  tau.I.save <- vector(length = n.mcmc)
  tau.P.save <- vector(length = n.mcmc)
  beta.0.save <- vector(length = n.mcmc)
  beta.1.save <- vector(length = n.mcmc)
  sigma.squared.save <- vector(length = n.mcmc)
  
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0) cat(' ', l)
    
    ##
    ## sample X
    ##
    
    ## for t = 0
    A.chol <- chol(alpha^2 * Sigma.epsilon.inv + 1 / sigma.squared.0.tilde * I)
    b <- alpha * Sigma.epsilon.inv %*% (X[0 + 1,] - (1 - alpha) * mu * J) + mu.0.tilde / sigma.squared.0.tilde * J
    X[0 + 1, ] <- rMVN(A.chol, b)
    
    # for t = 1,...,T - 1
    for(i in 1:(t - 1)){ # note the offset since t = 0 is the first element of the matrix X
      A.chol <- chol(t(Ht[[i]]) %*% Sigma.inv[[i]] %*% Ht[[i]] + (alpha^2 + 1) * Sigma.epsilon.inv)
      b <- t(Ht[[i]]) %*% Sigma.inv[[i]] %*% (Wt[[i]] - Bt[[i]]) + Sigma.epsilon.inv %*% (alpha * (X[i + 2, ] + X[i, ]) + (1 - alpha)^2 * mu)
      X[i + 1, ] <- rMVN(A.chol, b)
    }
    
    # for t = T
    A.chol <- chol(t(Ht[[t]]) %*% Sigma.inv[[t]] %*% Ht[[t]] + Sigma.epsilon.inv)
    b <- t(Ht[[t]]) %*% Sigma.inv[[t]] %*% (Wt[[t]] - Bt[[t]]) + Sigma.epsilon.inv %*% (alpha * X[t, ] + (1 - alpha) * mu)
    X[t + 1, ] <- rMVN(A.chol, b)
      
    ##
    ## sample beta_0
    ##
    
    A.chol <-sqrt(MP / tau.squared.P + 1 / sigma.squared.beta.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- sum(1 / tau.squared.P * IPt[[i]] %*% (WPt[[i]] - beta.1 * HPt[[i]] %*% X[i + 1, ]))
    }
    
    b <- sum(tmp) + mu.beta.0 / sigma.squared.beta.0
    beta.0 <- rMVN(A.chol, b)
    Bt <- lapply(1:t, make.Bt, beta.0 = beta.0)
    
    ##
    ## sample beta_1
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(HPt[[i]] %*% X[i, ]) %*% (HPt[[i]] %*% X[i, ]) / tau.squared.P
    }
    A.chol <-sqrt(sum(tmp) + 1 / sigma.squared.beta.1)       
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(HPt[[i]] %*% X[i, ]) %*% (WPt[[i]] - beta.0) / tau.squared.P
    }
    b <- sum(tmp) + mu.beta.1 / sigma.squared.beta.1
    beta.1 <- rMVN(A.chol, b)
    Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
    
    ##
    ## sample mu
    ##
    
    A.chol <- sqrt((1 - alpha)^2 * (MP / tau.squared.P + MI / tau.squared.I) + 1 / sigma.squared.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- sum((1 - alpha)^2 * Sigma.epsilon.inv %*% (X[i + 1, ] - alpha * X[i, ]))
    }
    b <- sum(tmp) + mu.0 / sigma.squared.0
    mu <- rMVN(A.chol, b)
    
    ##
    ## sample alpha
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i, ] - mu) %*% Sigma.epsilon.inv %*% (X[i, ] - mu)
    }
    A.inv <- 1 / sum(tmp)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i, ] - mu) %*% Sigma.epsilon.inv %*% (X[i + 1, ] - mu)
    }
    b <- sum(tmp)
    ##
    ## work on updating this truncated normal sampler to be more efficient
    ##
#     alpha.star <- rMVN(A.chol, b)
#     while(alpha.star < 0 || alpha.star > 1){
#       alpha.star <- rMVN(A.chol, b)
#     }
#     alpha <- alpha.star
    alpha <- rtruncnorm(1, a = 0, b = 1, mean = A.inv * b, sd = sqrt(A.inv))
    
    ##
    ## sample tau.squared.I
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      if(NIt[i] != 0){
        tmp[i] <- t(WIt[[i]] - HIt[[i]] %*% X[i + 1, ]) %*% (WIt[[i]] - HIt[[i]] %*% X[i + 1, ])
      } else {
        tmp[i] <- 0
      }
    }
  tau.squared.I <- 1 / rgamma(1, alpha.I + MI / 2, beta.I + sum(tmp) / 2)
  for(i in 1:t){
    Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]),rep(tau.squared.P, NPt[i])))
    Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]),rep(1 / tau.squared.P, NPt[i])))
  }    

    ##
    ## sample tau.squared.P
    ##

    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(WPt[[i]] - HPt[[i]] %*% X[i + 1, ] - beta.0) %*% (WPt[[i]] - HPt[[i]] %*% X[i + 1, ] - beta.0)
    }
    tau.squared.P <- 1 / rgamma(1, alpha.P + MP / 2, beta.I + sum(tmp) / 2)
    for(i in 1:t){
      Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]),rep(tau.squared.P, NPt[i])))
      Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]),rep(1 / tau.squared.P, NPt[i])))
    } 

    ##
    ## sample sigma.squared
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu) %*% Sigma.epsilon.inv %*% (X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu)
    }
    sigma.squared <- 1 / rgamma(1, alpha.sigma.squared + n * t / 2, beta.sigma.squared + sum(tmp) / 2)
    Sigma.epsilon <- sigma.squared * I
    Sigma.epsilon.inv <- 1 / sigma.squared * I

    ##
    ## save variables
    ##
    
    X.save[, , l] <- X
    alpha.save[l] <- alpha
    mu.save[l] <- mu
    phi.save[l] <- phi
    tau.I.save[l] <- tau.squared.I
    tau.P.save[l] <- tau.squared.P
    beta.0.save[l] <- beta.0
    beta.1.save[l] <- beta.1
    sigma.squared.save[l] <- sigma.squared
  }    
  
  list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save)  
}
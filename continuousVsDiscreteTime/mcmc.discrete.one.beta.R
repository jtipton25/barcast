##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

mcmc.disc <- function(WI, WP, HI, HP, params){
  
  ##
  ## libraries and subroutines
  ##
  
  make.neighborhood <- function(num.neighbors, dimension){ ## creates a neighborhood matrix for iCAR
    tmp <- matrix(0, dimension, dimension)
    for(i in 1:dimension){
      for(j in 1:dimension){
        if(i == j){
          tmp[i, j] <- 0
        } else if(abs(i - j) <= num.neighbors){
          tmp[i, j] <- 1
        }
      }
    }
    return(tmp)
  }
  
  make.tau.squared.P <- function(i, beta.0, beta.1, T){
    tmp <- (WP[i, ] - beta.1 * HP[i, ] * T[i] - beta.0 * HP[i, ])
    return(t(tmp) %*% tmp)
  }
  
  ##
  ## initialize variables
  ##
  
  n.mcmc <- params$n.mcmc
  delta.0 <- params$delta.0
  detla.1 <- params$delta.1
  num.truncate <- params$num.truncate
  alpha.I <- params$alpha.I
  beta.I <- params$beta.I
  alpha.P <- params$alpha.P
  beta.P <- params$beta.P
  alpha.sigma <- params$alpha.sigma
  beta.sigma <- params$beta.sigma
  num.neighbors <- params$num.neighbors
  
  t <- length(WI)
  tHP <- t(HP)
  H <- cbind(HI, HP)
  
  p <- dim(WP)[2]
  n.o <- length(which(WI != 0))
  n.u <- length(which(WI == 0))
  W <- cbind(WI, WP)
  Delta.0 <- diag(delta.0)
  Delta.1 <- diag(delta.1)
  
  NI.t <- vector(length = t)
  NP.t <- vector(length = t)
  for(i in 1:t){
    NI.t[i] <- HI[i]
    NP.t[i] <- sum(HP[i, ])
  }
  NI <- sum(NI.t)  
  NP <- sum(NP.t)
  
  ## initialize covariance matrix
  A <- make.neighborhood(num.neighbors, t)
  D <- diag(rowSums(A))
  Q <- D - A ## ICAR covariance matrix
  pc <- prcomp(Q)
  Z <- pc$rotation[, 1:num.truncate]
  tZ <- t(Z)
  tZZ <- tZ %*% Z
  Lambda <- diag((pc$sdev^2)[1:num.truncate])
  Lambda.inv <- diag(1 / (pc$sdev^2)[1:num.truncate])
  
  ## Initial latent field
  alpha <- rMVN(chol(Lambda.inv), rep(0, num.truncate))
  T <- Z %*% alpha
  T[HI == 1] <- WI[HI == 1]  ## initialize the latent field with the observed measurements
  
  ## initialize variance values 
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
  tmp <- T - Z %*% alpha	
  sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
  rm(tmp)
  
  ## initialize regression parameters
  beta.0 <- 0
  #   beta.0 <- rMVN(chol(tau.squared.P * solve(Delta.0)), rep(0, p))
  #   beta.0 <- rep(0, p)
  beta.1 <- rMVN(chol(tau.squared.P * 1 / Delta.1[1]), 0)
  J <- rep(1, p)
  #   beta.1 <- rMVN(chol(tau.squared.P * solve(Delta.1)), rep(0, p))
  
  ##
  ## set up save variables
  ##
  
  T.save <- matrix(0, nrow = t, ncol = n.mcmc)
  tau.squared.I.save <- rep(0, n.mcmc)
  tau.squared.P.save <- rep(0, n.mcmc)
  beta.0.save <- rep(0, n.mcmc)
  #   beta.0.save <- matrix(0, nrow = p, ncol = n.mcmc)
  beta.1.save <- rep(0, n.mcmc)
  #   beta.1.save <- matrix(0, nrow = p, ncol = n.mcmc)
  sigma.squared.save <- rep(0, n.mcmc)
  trend.save <- matrix(0, nrow = t, ncol = n.mcmc)
  alpha.save <- matrix(0, nrow = num.truncate, ncol = n.mcmc)
  
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0){
      cat(' ', l)
    } 
    
    ##
    ## sample latent field T
    ##
    
    b <- rep(0, t)
    A.tmp.vec <- rep(0, t)
    for(i in 1:t){
      A.tmp.vec[i] <- (H[i, ] * c(1, beta.1 * J)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1 * J))
      b[i] <- (H[i, ] * c(1, beta.1 * J)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0 * J))
      #     	A.tmp.vec[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1))
      #       b[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0))
    }
    A.tmp.chol <- chol(diag(A.tmp.vec + 1 / sigma.squared))
    T <- rMVN(A.tmp.chol, b + Z %*% alpha / sigma.squared)
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample beta_0
    ##
    
    A.tmp <- 0
    b <- 0
    #     A.tmp <- matrix(0, p, p)
    #     b <- rep(0, p)
    for(i in 1:t){
      if(HI[i] == 1){
        A.tmp <- A.tmp + HP[i, ] %*% HP[i, ]
        #         A.tmp <- A.tmp + HP[i, ] %*% t(HP[i, ])
        b <- b + HP[i, ] %*% (WP[i, ] - beta.1 * HP[i, ] * T[i])
      }
    }
    beta.0 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.0[1])), b / tau.squared.P)
    #     beta.0 <- rMVN(chol(1 / tau.squared.P * A.tmp + Delta.0), b / tau.squared.P)
    rm(A.tmp)    
    rm(b)
    
    ##
    ## sample beta_1
    ##
    
    A.tmp <- 0
    b <- 0
    #     A.tmp <- matrix(0, p, p)
    #     b <- rep(0, p)
    for(i in 1:t){ 
      if(HI[i] == 1){ ## only sample for observed data???
        A.tmp <- A.tmp + T[i]^2 * HP[i, ] %*% HP[i, ]
        #        A.tmp <- A.tmp + T[i]^2 * HP[i, ] %*% t(HP[i, ])
        b <- b + T[i] * HP[i, ] %*% (WP[i, ] - beta.0 * HP[i, ])
      }
    }
    beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1[1])), b / tau.squared.P)
    #     beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1)), b / tau.squared.P)
    rm(A.tmp)
    rm(b)    
    
    ##
    ## sample tau.squared.I
    ##
    
    tau.squared.I <- 1 / rgamma(1, alpha.I + NI / 2, beta.I + sum((WI * HI - T * HI)^2) / 2)
    
    ##
    ## sample tau.squared.P
    ##
    
    tau.squared.P <- 1 / rgamma(1, alpha.P + NP / 2, beta.P + sum(sapply(1:t, make.tau.squared.P, beta.0 = beta.0, beta.1 = beta.1, T = T)) / 2)
    Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
    
    ##
    ## sample alpha
    ##
    
    A.tmp.chol <- chol(tZZ / sigma.squared + Lambda.inv)
    b <- tZ %*% T / sigma.squared
    alpha <- rMVN(A.tmp.chol, b)	
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- T - Z %*% alpha	
    sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
    rm(tmp)
    
    ##
    ## save variables
    ##
    
    T.save[, l] <- T
    tau.squared.I.save[l] <- tau.squared.I
    tau.squared.P.save[l] <- tau.squared.P
    beta.0.save[l] <- beta.0
    #     beta.0.save[, l] <- beta.0
    beta.1.save[l] <- beta.1
    #     beta.1.save[, l] <- beta.1
    sigma.squared.save[l] <- sigma.squared
    alpha.save[, l] <- alpha
    trend.save[, l]<- Z %*% alpha
  }
  list(T.save = T.save, tau.squared.I.save = tau.squared.I.save, tau.squared.P.save = tau.squared.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, alpha.save = alpha.save, trend.save = trend.save)
}
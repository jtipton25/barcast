##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

## Next: Remove Sigma, Sigma.inv, Sigma.epsilon and Sigma.epsilon.inv
##
##
## Note that in the code X represents the field T in the Barcast Model
mcmc <- function(WI, WP, HI, HP, params, method = 'continuous'){
  
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
  
  make.tau.squared.P <- function(i, beta.1, T){
    tmp <- (WP[i, ] - beta.1 * HP[i, ] * T[i])
    return(t(tmp) %*% tmp)
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
  if(method == 'discrete'){
    num.neighbors <- params$num.neighbors
  }
  if(method == 'continuous'){
    phi.squared.prior <- params$phi.squared.prior
  }
  
  t <- length(WI)
  tHP <- t(HP)
  H <- cbind(HI, HP)
  
  p <- dim(WP)[2]
  n.o <- length(which(WI != 0))
  n.u <- length(which(WI == 0))
  W <- cbind(WI, WP)
  Delta.0 <- diag(delta.0)
  Delta.1 <- diag(delta.1)
  
  ## initialize variance values <- this is naive, could do better?
  
  tau.squared.I <- 1 / rgamma(1, 1, 1)
  tau.squared.P <- 1 / rgamma(1, 1, 1)
  sigma.squared <- 1 / rgamma(1, 1, 1)
  Sigma <- diag(c(tau.squared.I, rep(tau.squared.P, p)))
  Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
  
  NI.t <- vector(length = t)
  NP.t <- vector(length = t)
  for(i in 1:t){
    NI.t[i] <- HI[i]
    NP.t[i] <- sum(HP[i, ])
  }
  NI <- sum(NI.t)  
  NP <- sum(NP.t)
  
  ## initialize covariance matrix
  if(method == 'continuous'){
    D <- as.matrix(dist(1:t))
    Q <- exp( - D^2 / phi.squared.prior) ## gaussian covariance matrix
  } else if(method == 'discrete'){
    A <- make.neighborhood(num.neighbors, t)
    D <- diag(rowSums(A))
    Q <- D - A ## ICAR covariance matrix
  } else {
    stop("not a valid method")
  }
  
  pc <- prcomp(Q)
  Z <- pc$rotation
  tZ <- t(Z)
  tZZ <- tZ %*% Z
  Lambda <- diag(pc$sdev^2)
  Lambda.inv <- diag(1 / pc$sdev^2)
  alpha <- rMVN(chol(Lambda.inv), rep(0, t))
  devs <- rnorm(t, 0, sigma.squared)
  
  T <- Z %*% alpha + devs
  T[HI == 1] <- WI[HI == 1] ## initialize the latent field with the observed measurements
  
  ## initialize regression parameters
  beta.0 <- rep(0, p)
  beta.1 <- rnorm(p)
  
  ##
  ## set up save variables
  ##
  
  T.save <- matrix(nrow = t, ncol = n.mcmc)
  tau.squared.I.save <- vector(length = n.mcmc)
  tau.squared.P.save <- vector(length = n.mcmc)
  beta.0.save <- matrix(nrow = p, ncol = n.mcmc)
  beta.1.save <- matrix(nrow = p, ncol = n.mcmc)
  sigma.squared.save <- vector(length = n.mcmc)
  alpha.save <- matrix(nrow = t, ncol = n.mcmc)
  trend.save <- matrix(nrow = t, ncol = n.mcmc)
    
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0) cat(' ', l)
    
    ##
    ## sample latent field T
    ##
    
    A.tmp.chol <- chol(diag(NI.t * 1 / tau.squared.I + NP.t * 1 / tau.squared.P + sigma.squared))
    b <- rep(0, t)
    for(i in 1:t){
      b[i] <- H[i, ] %*% Sigma.inv %*% W[i, ]
    }
    T <- rMVN(A.tmp.chol, b + Z %*% alpha / sigma.squared)
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample beta_0
    ##
    
    ## maybe we need to add a regularization... 
    A.tmp <- matrix(0, p, p)
    b <- rep(0, p)
    for(i in 1:t){
      if(HI[i] == 1){
        #         A.tmp <- A.tmp + HP[i, ] %*% tHP[, i]
        A.tmp <- A.tmp + HP[i, ] %*% t(HP[i, ])
        b <- b + HP[i, ] %*% (WP[i, ] - beta.1 * HP[i, ] * T[i])
      }
    }
    beta.0 <- rMVN(chol(1 / tau.squared.P * A.tmp + Delta.0), b / tau.squared.P)
    rm(A.tmp)    
    rm(b)
    
    ##
    ## sample beta_1
    ##
    
    A.tmp <- matrix(0, p, p)
    b <- rep(0, p)
    for(i in 1:t){
      if(HI[i] == 1){
        A.tmp <- A.tmp + T[i]^2 * HP[i, ] %*% t(HP[i, ])
        b <- b + T[i] * HP[i, ] * (WP[i, ] - beta.0 * HP[i, ])
      }
    }
    beta.1 <- rMVN(chol(1 / tau.squared.P * A.tmp + Delta.1), b / tau.squared.P)
    rm(A.tmp)
    rm(b)
    
    
    ##
    ## sample tau.squared.I
    ##
    
    tau.squared.I <- 1 / rgamma(1, NI / 2, sum((WI * HI - T * HI)^2) / 2)
    
    ##
    ## sample tau.squared.P
    ##
    
    tau.squared.P <- 1 / rgamma(1, NP / 2, sum(sapply(1:t, make.tau.squared.P, beta.0 = beta.0, beta.1 = beta.1, T = T)) / 2)
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
    sigma.squared <- 1 / rgamma(1, t / 2, t(tmp) %*% tmp / 2)
    rm(tmp)
    
    ##
    ## save variables
    ##
    
    T.save[, l] <- T
    tau.squared.I.save[l] <- tau.squared.I
    tau.squared.P.save[l] <- tau.squared.P
    beta.0.save[, l] <- beta.0
    beta.1.save[, l] <- beta.1
    sigma.squared.save[l] <- sigma.squared
    alpha.save[, l] <- alpha
    trend.save[, l] <- Z %*% alpha
  }    
  list(T.save = T.save, tau.squared.I.save = tau.squared.I.save, tau.squared.P.save = tau.squared.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, alpha.save = alpha.save, trend.save = trend.save)
}
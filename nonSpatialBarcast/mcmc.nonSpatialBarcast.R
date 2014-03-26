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
  
  ##
  ## initialize variables
  ##
  
  W <- as.matrix(cbind(WI, WP))
  WP <- as.matrix(WP)
  t <- dim(W)[1]
  n <- dim(W)[2]
#   t.u <- c(1:451, 556) #years PDSI is unobserved
#   t.o <- 452:555 # years PDSI is observed
  HI <- !is.na(WI)
  HP <- !is.na(WP)
  H <- !is.na(W) 
  X <- matrix(0, nrow = t + 1, ncol = n)
  NI <- 1
  NP <- dim(WP)[2]
  NIt <- vector(length = t)
  NPt <- vector(length = t)
  for(i in 1:t){
    NIt[i] <- sum(!is.na(WI[i]))
    NPt[i] <- sum(!is.na(WP[i, ]))
  }
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

  Sigma.epsilon <- sigma.squared * I
  Sigma.epsilon.inv <- 1 / sigma.squared * I
#   Sigma <- vector('list', length = t)
#   Sigma.inv <- vector('list', length = t)
#   for(i in 1:t){
#     Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]),rep(tau.squared.P, NPt[i])))
#     Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]),rep(1 / tau.squared.P, NPt[i])))
#   }  
  Sigma <- c(rep(tau.squared.I, NI), rep(tau.squared.P, NP)) * I
  Sigma.inv <- c(rep(1 / tau.squared.I, NI), rep(1 / tau.squared.P, NP)) * I   
  B <- c(rep(0, NI), rep(beta.0, NP))

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
  
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0) cat(' ', l)
    
    ##
    ## sample X
    ##
    
    ## for t = 0
    A.chol <- chol(alpha * Sigma.epsilon.inv + 1 / sigma.squared.0.tilde * I)
    b <- alpha * Sigma.inv %*% X[1,] - (1 - alpha) * mu + mu.0.tilde / sigma.squared.0.tilde  
    X[1, ] <- rMVN(A.chol, b)
    
    # for t = 1,...,T - 1
    for(i in 2:t){ # note the offset since t = 0 is the first element of the matrix X
      A.chol <- chol(Sigma.inv + (alpha^2 + 1) * Sigma.epsilon.inv)
      b <- Sigma.inv[, H[i - 1, ]] %*% (W[i - 1, ][H[i - 1, ]] - B[H[i - 1, ]]) + Sigma.epsilon.inv %*% (alpha * (X[i + 1, ] + X[i - 1, ]) + (1 - alpha)^2 * mu)
      X[i, ] <- rMVN(A.chol, b)
    }
    
    # for t = T
    A.chol <- chol(Sigma.inv + Sigma.epsilon.inv)
    b <- Sigma.inv[, H[t, ]] %*% (W[t, ][H[t, ]] - B[H[t, ]]) + Sigma.epsilon.inv %*% (alpha * X[t - 1, ] + (1 - alpha) * mu)
    X[t + 1, ] <- rMVN(A.chol, b)
      
    ##
    ## sample beta_0
    ##
    
    A.chol <-sqrt(MP / tau.squared.P + 1 / sigma.squared.beta.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- sum(1 / tau.squared.P * IPt[[i]] %*% (WP[i, ][HP[i, ]] - beta.1 * X[i, 2:n][HP[i, ]]))
    }
    b <- sum(tmp) + mu.beta.0 / sigma.squared.beta.0
    beta.0 <- rMVN(A.chol, b)
    
    ##
    ## sample beta_1
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i, 2:n][HP[i, ]]) %*% X[i, 2:n][HP[i, ]] / tau.squared.P
    }
    A.chol <-sqrt(sum(tmp) + 1 / sigma.squared.beta.1)       
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i, 2:n][HP[i, ]]) %*% (1 / tau.squared.P * IPt[[i]] %*% (WP[i, ][HP[i, ]] - beta.0))
    }
    b <- sum(tmp) + mu.beta.1 / sigma.squared.beta.1
    beta.1 <- rMVN(A.chol, b)
    
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
    
    tau.squared.I <- 1 / rgamma(1, alpha.I + MI / 2, beta.I + t(WI[HI] - X[, 1][HI]) %*% (WI[HI] - X[, 1][HI]) / 2)

    ##
    ## sample tau.squared.P
    ##

    tau.squared.P <- 1 / rgamma(1, alpha.P + MP / 2, beta.I + t(WP[HP] - X[, 2:n][HP]) %*% (WP[HP] - X[, 2:n][HP]) / 2)
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu) %*% Sigma.epsilon.inv %*% (X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu)
    }
    sigma.squared <- 1 / rgamma(1, alpha.sigma.squared + n * t / 2, beta.sigma.squared + sum(tmp) / 2)
    
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
  }    
  
  list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save)  
}
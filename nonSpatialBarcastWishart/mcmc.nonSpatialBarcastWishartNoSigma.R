##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

## Next: Remove Sigma, Sigma.inv, Sigma.epsilon and Sigma.epsilon.inv
##
##
## Note that in the code X represents the field T in the Barcast Model
mcmc <- function(WI, WP, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.sigma.squared, beta.sigma.squared, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1, nu.wish){
  
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
  
  make.Wt <- function(i){
    as.vector(W[i, ][H[i,]])
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

  make.mh <- function(i, Z, mu, alpha, Sigma.epsilon, Sigma.epsilon.inv){
    - 1 / 2 * determinant(Sigma.epsilon, logarithm = TRUE)$modulus[1] - 1 / 2 * t(Z[i + 1, ] - alpha * Z[i, ] - (1 - alpha) * mu) %*% Sigma.epsilon.inv %*% (Z[i + 1, ] - alpha * Z[i, ] - (1 - alpha) * mu)
  }
  
  ##
  ## initialize variables
  ##
  
  W <- as.matrix(cbind(WI, WP))
  WP <- as.matrix(WP)
  t <- dim(W)[1]
  n <- dim(W)[2]
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
  MI <- sum(NIt)  
  MP <- sum(NPt)
  
  HI.idx <- !is.na(WI)
  HI <- ifelse(HI.idx, 1, 0)
  HP.idx <- !is.na(WP)
  HP <- ifelse(HP.idx, 1, 0)
  H <- !is.na(W)
  
  Wt <- lapply(1:t, make.Wt)
  WPt <- lapply(1:t, make.WPt)
  WIt <- lapply(1:t, make.WIt) 
  HIt <- lapply(1:t, make.HIt)
  HPt <- lapply(1:t, make.HPt)  
  J <- rep(1, n)
  I <- diag(n)
  
  
  alpha <- runif(1, alpha.alpha, beta.alpha)
  mu <- rnorm(1, mu.0, sigma.squared.0)
#   sigma.squared <- 1 / rgamma(1, alpha.sigma.squared, beta.sigma.squared)
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  beta.1 <- rnorm(1, mu.beta.1, sigma.squared.beta.1)
  Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
  
  Q <- riwish(nu.wish, I) # * (nu.wish - n - 1)
  Q.inv <- solve(Q)
#   Sigma.epsilon <- sigma.squared * Q
  Sigma.epsilon <- Q
#   Sigma.epsilon.inv <- 1 / sigma.squared * Q.inv
  Sigma.epsilon.inv <- Q.inv
  Sigma.full <- diag(c(rep(tau.squared.I, NI), rep(tau.squared.P, NP)))
  Sigma.full.inv <- diag(c(rep(1 / tau.squared.I, NI), rep(1 / tau.squared.P, NP)))
  Sigma <- vector('list', length = t)
  Sigma.inv <- vector('list', length = t)
  for(i in 1:t){
    Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]), rep(tau.squared.P, NPt[i])))
    Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]), rep(1 / tau.squared.P, NPt[i])))
  }  
  
  ##
  ## set up save variables
  ##
  
  X.save <- array(dim = c(t + 1, n, n.mcmc))
  alpha.save <- vector(length = n.mcmc)
  mu.save <- vector(length = n.mcmc)
  phi.save <- vector(length = n.mcmc)
  tau.I.save <- vector(length = n.mcmc)
  tau.P.save <- vector(length = n.mcmc)
  beta.1.save <- vector(length = n.mcmc)
#   sigma.squared.save <- vector(length = n.mcmc)
  Q.save <- array(dim = c(n, n, n.mcmc))
  
  ##
  ## mcmc loop
  ##
  
  for(l in 1:n.mcmc){
    if(l %% 100 == 0) cat(' ', l)
    
    ##
    ## sample X
    ##
    
    ## for t = 0
    A.chol <- chol(alpha^2 * Sigma.epsilon.inv + I / sigma.squared.0.tilde)
    b <- alpha * Sigma.epsilon.inv %*% (X[2, ] - (1 - alpha) * mu) + mu.0.tilde / sigma.squared.0.tilde
    X[1, ] <- rMVN(A.chol, b)
    
    # for t = 1,...,T - 1
    for(i in 1:(t - 1)){ # note the offset since t = 0 is the first element of the matrix X
      A.chol <- chol(t(Ht[[i]]) %*%  Sigma.inv[[i]] %*% Ht[[i]] + (alpha^2 + 1) * Sigma.epsilon.inv)
      b <- t(Ht[[i]]) %*% Sigma.inv[[i]] %*% (Wt[[i]]) + Sigma.epsilon.inv %*% (alpha * (X[i + 2, ] + X[i, ]) + (1 - alpha)^2 * mu)
      X[i + 1, ] <- rMVN(A.chol, b)
    }
    
    # for t = T
    A.chol <- chol(t(Ht[[t]]) %*% Sigma.inv[[t]] %*% Ht[[t]] + Sigma.epsilon.inv)
    b <- t(Ht[[t]]) %*% Sigma.inv[[t]] %*% (Wt[[t]]) + Sigma.epsilon.inv %*% (alpha * X[t, ] + (1 - alpha) * mu)
    X[t + 1, ] <- rMVN(A.chol, b)
    
    ##
    ## sample beta_1
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(HPt[[i]] %*% X[i + 1, ]) %*% (HPt[[i]] %*% X[i + 1, ])
    }
    A.chol <-sqrt(sum(tmp) / tau.squared.P + 1 / sigma.squared.beta.1)       
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- t(HPt[[i]] %*% X[i + 1, ]) %*% (WPt[[i]])
    }
    b <- sum(tmp) / tau.squared.P + mu.beta.1 / sigma.squared.beta.1
    beta.1 <- rMVN(A.chol, b)
    Ht <- lapply(1:t, make.Ht, beta.1 = beta.1)
    
    ##
    ## sample mu
    ##
    
    #     A.chol <- sqrt((1 - alpha)^2 * (MP / tau.squared.P + MI / tau.squared.I) + 1 / sigma.squared.0)
    A.chol <- sqrt(t * (1 - alpha)^2 * J %*% Sigma.epsilon.inv %*% J + 1 / sigma.squared.0)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- J %*% Sigma.epsilon.inv %*% (X[i + 1, ] - alpha * X[i, ])
    }
    b <- (1 - alpha) * sum(tmp) + mu.0 / sigma.squared.0
    mu <- rMVN(A.chol, b)
    
    ##
    ## sample alpha
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- (X[i, ] - mu) %*% Sigma.epsilon.inv %*% (X[i, ] - mu)
    }
    A.inv <- 1 / sum(tmp)
    tmp <- vector(length = t)
    for(i in 1:t){
      tmp[i] <- (X[i, ] - mu) %*% Sigma.epsilon.inv %*% (X[i + 1, ] - mu)
    }
    b <- sum(tmp)
    alpha <- rtruncnorm(1, a = 0, b = 1, mean = A.inv * b, sd = sqrt(A.inv))
    
    ##
    ## sample tau.squared.I
    ##
    
    tmp <- vector(length = t)
    for(i in 1:t){
      if(NIt[i] != 0){
        tmp[i] <- (WIt[[i]] - HIt[[i]] %*% X[i + 1, ]) %*% (WIt[[i]] - HIt[[i]] %*% X[i + 1, ])
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

#       tmp[i] <- t(WPt[[i]] - beta.1 * HPt[[i]] %*% X[i + 1, ] - beta.0) %*% (WPt[[i]] - beta.1 * HPt[[i]] %*% X[i + 1, ] - beta.0)
      tmp[i] <- t(WPt[[i]] - beta.1 * HPt[[i]] %*% X[i + 1, ]) %*% (WPt[[i]] - beta.1 * HPt[[i]] %*% X[i + 1, ])
    }
    tau.squared.P <- 1 / rgamma(1, alpha.P + MP / 2, beta.P + sum(tmp) / 2)
    
    for(i in 1:t){
      Sigma[[i]] <- diag(c(rep(tau.squared.I, NIt[i]),rep(tau.squared.P, NPt[i])))
      Sigma.inv[[i]] <- diag(c(rep(1 / tau.squared.I, NIt[i]),rep(1 / tau.squared.P, NPt[i])))
    } 
    
    ##
    ## sample sigma.squared
    ##
#     
#     tmp <- vector(length = t)
#     for(i in 1:t){
#       tmp[i] <- (X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu) %*% Q.inv %*% (X[i + 1, ] - alpha * X[i, ] - (1 - alpha) * mu)
#     }
#     sigma.squared <- 1 / rgamma(1, alpha.sigma.squared + n * t / 2, beta.sigma.squared + sum(tmp) / 2)
#     Sigma.epsilon <- sigma.squared * Q
#     Sigma.epsilon.inv <- 1 / sigma.squared * Q.inv
    
    ## 
    ## sample Q
    ##
    
    
#     Tbar <- apply(X[ - 1, ], 1, mean) - alpha * apply(X[ - (T + 1), ], 1, mean)
    Tbar <- apply(X[ - 1, ], 2, mean) - alpha * apply(X[ - (T + 1), ], 2, mean)
    Tbar.mat <- matrix(rep(Tbar, t), nrow = t, ncol = n, byrow = TRUE)
    QT <- (X[ - 1, ] -  alpha * X[ - (t + 1), ] - Tbar.mat)
    Q0 <- t(QT) %*% QT
    q <- (Tbar - (1 - alpha) * mu)
#     Q <- riwish(t + nu.wish, (Q0 + t * q %*% t(q) + I) / sigma.squared)
    Q <- riwish(t + nu.wish, (Q0 + t * q %*% t(q) + I)) # * (nu.wish + t - n - 1)
    Q.inv <- solve(Q)
#     Sigma.epsilon. <- sigma.squared * Q
#     Sigma.epsilon.inv <- 1 / sigma.squared * Q.inv
    Sigma.epsilon. <- Q
    Sigma.epsilon.inv <- Q.inv

    ##
    ## save variables
    ##
    
    X.save[, , l] <- X
    alpha.save[l] <- alpha
    mu.save[l] <- mu
    tau.I.save[l] <- tau.squared.I
    tau.P.save[l] <- tau.squared.P
    beta.1.save[l] <- beta.1
#     sigma.squared.save[l] <- sigma.squared
    Q.save[, , l] <- Q
  }    
  
#   list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, rho.save = rho.save, rho.accept = rho.accept)  
#   list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, Q.save = Q.save)  
  list(X.save = X.save, alpha.save = alpha.save, mu.save = mu.save, phi.save = phi.save, tau.I.save = tau.I.save, tau.P.save = tau.P.save, beta.1.save = beta.1.save, Q.save = Q.save)  
}
##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

mcmc.cont <- function(WI, WP, HI, HP, params){
  
  ##
  ## libraries and subroutines
  ##
  
  make.tau.squared.P <- function(i, beta.0, beta.1, T){
    ## only calculate for observed T
#     if(HI[i] == 0){
#       tmp <- 0
#     } else {
      tmp <- (WP[i, ] - beta.1 * HP[i, ] * T[i] - beta.0 * HP[i, ])
#     }
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
  phi.lower <- params$phi.lower
  phi.upper <- params$phi.upper
  N.phi <- params$N.phi
  phi.tune <- params$phi.tune
  phi.prior <- seq(phi.lower, phi.upper, length = N.phi)
  
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
  D <- as.matrix(dist(1:t))
  Z <- vector('list', length = N.phi)
  tZ <- vector('list', length = N.phi)
  Lambda <- vector('list', length = N.phi)
  Lambda.inv <- vector('list', length = N.phi)
  for(i in 1:N.phi){
    Q <- exp( - D / phi.prior[i]) ## gaussian covariance matrix
    pc <- prcomp(Q)
#     Z[[i]] <- pc$rotation[, 1:num.truncate]
		Z[[i]] <- pc$rotation
    tZ[[i]] <- t(Z[[i]])
#     Lambda[[i]] <- diag((pc$sdev^2)[1:num.truncate])
#     Lambda.inv[[i]] <- diag((1 / pc$sdev^2)[1:num.truncate])
    Lambda[[i]] <- diag(pc$sdev^2)
    Lambda.inv[[i]] <- diag(1 / pc$sdev^2)
  }
  
  I.t <- diag(rep(1, t))
  I.trunc <- diag(rep(1, num.truncate))
  phi.idx <- sample(1:N.phi, 1)
  
  ## initialize latent field
  alpha <- rMVN(chol(Lambda.inv[[phi.idx]]), rep(0, num.truncate))
  T <- Z[[i]] %*% alpha
  T[HI == 1] <- WI[HI == 1]  ## initialize the latent field with the observed measurements
  
  ## initialize variance values 
#   tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
tau.squared.I <- 0.1
#   tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
tau.squared.P <- 1.25
  Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
  tmp <- T - Z[[phi.idx]] %*% alpha	
#   sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
sigma.squared <- 0.5
  rm(tmp)
  
  ## initialize regression parameters
  #   beta.0 <- rMVN(chol(tau.squared.P * solve(Delta.0)), rep(0, p))
  beta.0 <- 0
  #   beta.1 <- rMVN(chol(tau.squared.P * solve(Delta.1)), rep(0, p))
  beta.1 <- rMVN(chol(tau.squared.P * 1 / Delta.1), 0)
	J <- rep(1, p)
	beta.mat <- matrix(c(1, beta.1 * J), nrow = t, ncol = p + 1, byrow = TRUE)
  
  
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
  phi.idx.save <- rep(0, n.mcmc)
  phi.accept <- 0
  
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
    
    
    
#     b <- rep(0, t)
#     A.tmp.vec <- rep(0, t)
#     for(i in 1:t){
#       A.tmp.vec[i] <- (H[i, ] * c(1, beta.1 * J)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1 * J))
#       b[i] <- (H[i, ] * c(1, beta.1 * J)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0 * J))
      #       A.tmp.vec[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1))
      #       b[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0)
#     }
#     A.tmp.chol <- chol(diag(A.tmp.vec + 1 / sigma.squared))
#     test <- 
#     chol((rowSums((H * beta.mat) %*% Sigma.inv * (H * beta.mat)) + 1 / sigma.squared) * I.t)
#     A.tmp.vec

    T <- rMVN(chol((rowSums((H * beta.mat) %*% Sigma.inv * (H * beta.mat)) + 1 / sigma.squared) * I.t), (rowSums((H * beta.mat) %*% Sigma.inv * (W)) + Z[[phi.idx]] %*% alpha / sigma.squared))
    T.mat <- matrix(T, nrow = t, ncol = p, byrow = FALSE)
#     rm(A.tmp.vec)
#     rm(A.tmp.chol)
#     rm(b)
    
    ##
    ## sample beta_0
    ##
    
#     A.tmp <- 0
#     b <- 0
#     #         A.tmp <- matrix(0, p, p)
#     #         b <- rep(0, p)
#     for(i in 1:t){
#       if(HI[i] == 1){ ## only sample for observed data???
#         A.tmp <- A.tmp + HP[i, ] %*% HP[i, ]
#         #             A.tmp <- A.tmp + HP[i, ] %*% t(HP[i, ])
#         b <- b + HP[i, ] %*% (WP[i, ] - beta.1 * HP[i, ] * T[i])
#       }
#     }
#     beta.0 <- rMVN(chol(1 / tau.squared.P *(A.tmp + Delta.0[1])), b / tau.squared.P)
#     rm(A.tmp)    
#     rm(b)
    
    ##
    ## sample beta_1
    ##
    
#     A.tmp <- 0
#     b <- 0
#     #     A.tmp <- matrix(0, p, p)
#     #     b <- rep(0, p)
#     for(i in 1:t){ 
#       if(HI[i] == 1){ ## only sample for observed data???
#         A.tmp <- A.tmp + T[i]^2 * HP[i, ] %*% HP[i, ]
#         #        A.tmp <- A.tmp + T[i]^2 * HP[i, ] %*% t(HP[i, ])
#         b <- b + T[i] * HP[i, ] %*% (WP[i, ] - beta.0 * HP[i, ])
#       }
#     }

# 		sum(NP.t[t.o] * T[t.o]^2)
# 		sum(NP * T^2)
#     beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1[1])), b / tau.squared.P)
    #     beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1)), b / tau.squared.P)
# 		beta.1 <- rMVN(chol((sum(NP.t[t.o] * T[t.o]^2) + Delta.1) / tau.squared.P), sum(apply((T.mat * HP) * (WP), 1, sum)[t.o]) / tau.squared.P)
		beta.1 <- rMVN(chol((sum(NP.t * T^2) + Delta.1) / tau.squared.P), sum((T.mat * HP) * (WP)) / tau.squared.P)	 
#     rm(A.tmp)
#     rm(b)
		beta.mat <- matrix(c(1, beta.1 * J), nrow = t, ncol = p + 1, byrow = TRUE)
    
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
    
    A.tmp.chol <- chol(I.trunc / sigma.squared + Lambda.inv[[phi.idx]])
    b <- tZ[[phi.idx]] %*% T / sigma.squared
    alpha <- rMVN(A.tmp.chol, b)
		talpha <- t(alpha)
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- T - Z[[phi.idx]] %*% alpha	
    sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
    rm(tmp)
    
    ##
    ## sample phi for continuous time model
    ##
    
		phi.idx.star <- phi.idx + sample( (- phi.tune:phi.tune)[ - (phi.tune + 1)], 1)
		if(phi.idx.star > N.phi){
			phi.idx.star <- N.phi
		}
		if(phi.idx.star < 1){
			phi.idx.star <- 1
		}
      mh1 <- - 1 / (2 * sigma.squared) * ( - 2 * t(T) %*% Z[[phi.idx.star]] %*% alpha + talpha %*% alpha) - 1 / 2 * talpha %*% Lambda.inv[[phi.idx.star]] %*% alpha
 			mh2 <- - 1 / (2 * sigma.squared) * ( - 2 * t(T) %*% Z[[phi.idx]] %*% alpha + talpha %*% alpha) - 1 / 2 * talpha %*% Lambda.inv[[phi.idx]] %*% alpha
#       rm(tmp)
      mh <- exp(mh1 - mh2)
      if(mh > runif(1)){
        phi.idx <- phi.idx.star
        phi.accept <- phi.accept + 1 / n.mcmc
      }
    
    ##
    ## save variables
    ##
    
    T.save[, l] <- T
    tau.squared.I.save[l] <- tau.squared.I
    tau.squared.P.save[l] <- tau.squared.P
    beta.0.save[l] <- beta.0
    #     beta.0.save[, l] <- beta.0
    #     beta.1.save[, l] <- beta.1
    beta.1.save[l] <- beta.1
    sigma.squared.save[l] <- sigma.squared
    alpha.save[, l] <- alpha
    trend.save[, l] <- Z[[phi.idx]] %*% alpha
    phi.idx.save[l] <- phi.idx
  }
  list(T.save = T.save, tau.squared.I.save = tau.squared.I.save, tau.squared.P.save = tau.squared.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, alpha.save = alpha.save, trend.save = trend.save, phi.accept = phi.accept, phi.idx.save = phi.idx.save)
}
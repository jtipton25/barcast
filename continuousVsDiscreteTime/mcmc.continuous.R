##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

mcmc.cont <- function(WI, WP, HI, HP, params){
  
  ##
  ## libraries and subroutines
  ##
  
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
  phi.squared.lower <- params$phi.squared.lower
  phi.squared.upper <- params$phi.squared.upper
  N.phi <- params$N.phi
  phi.tune <- params$phi.tune
  phi.squared.prior <- seq(phi.squared.lower, phi.squared.upper, length = N.phi)

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
  tZZ <- vector('list', length = N.phi)
  Lambda <- vector('list', length = N.phi)
  Lambda.inv <- vector('list', length = N.phi)
  for(i in 1:N.phi){
    Q <- exp( - D^2 / phi.squared.prior[i]) ## gaussian covariance matrix
    pc <- prcomp(Q)
    Z[[i]] <- pc$rotation[, 1:num.truncate]
    tZ[[i]] <- t(Z[[i]])
    tZZ[[i]] <- tZ[[i]] %*% Z[[i]]
    Lambda[[i]] <- diag((pc$sdev^2)[1:num.truncate])
    Lambda.inv[[i]] <- diag((1 / pc$sdev^2)[1:num.truncate])
  }
	
	phi.idx <- sample(1, 1:N.phi)
	alpha <- rMVN(chol(Lambda.inv[[phi.idx]]), rep(0, num.truncate))

	T <- Z[[i]] %*% alpha
  T[HI == 1] <- WI[HI == 1]  ## initialize the latent field with the observed measurements
    
  ## initialize variance values 
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
  tmp <- T - Z[[phi.idx]] %*% alpha	
  sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
  rm(tmp)
  
  ## initialize regression parameters
#   beta.0 <- rMVN(chol(tau.squared.P * solve(Delta.0)), rep(0, p))
  beta.0 <- rep(0, p)
  beta.1 <- rMVN(chol(tau.squared.P * solve(Delta.1)), rep(0, p))

    
  ##
  ## set up save variables
  ##
  
  T.save <- matrix(0, nrow = t, ncol = n.mcmc)
  tau.squared.I.save <- rep(0, n.mcmc)
  tau.squared.P.save <- rep(0, n.mcmc)
  beta.0.save <- matrix(0, nrow = p, ncol = n.mcmc)
  beta.1.save <- matrix(0, nrow = p, ncol = n.mcmc)
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
    
#     A.tmp.chol <- matrix(0, nrow = t, ncol = t)
#     	chol(diag(NI.t * 1 / tau.squared.I + NP.t * 1 / tau.squared.P + sigma.squared))
    b <- rep(0, t)
		A.tmp.vec <- rep(0, t)
    for(i in 1:t){
    	A.tmp.vec[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1))
      b[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0))
    }
		A.tmp.chol <- chol(diag(A.tmp.vec + 1 / sigma.squared))
    T <- rMVN(A.tmp.chol, b + Z[[phi.idx]] %*% alpha / sigma.squared)
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample beta_0
    ##
    
#     A.tmp <- matrix(0, p, p)
#     b <- rep(0, p)
#     for(i in 1:t){
#       if(HI[i] == 1){
#         #         A.tmp <- A.tmp + HP[i, ] %*% tHP[, i]
#         A.tmp <- A.tmp + HP[i, ] %*% t(HP[i, ])
#         b <- b + HP[i, ] %*% (WP[i, ] - beta.1 * HP[i, ] * T[i])
#       }
#     }
#     beta.0 <- rMVN(chol(1 / tau.squared.P * A.tmp + Delta.0), b / tau.squared.P)
#     rm(A.tmp)    
#     rm(b)
#     
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
    beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1)), b / tau.squared.P)
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
    
  	A.tmp.chol <- chol(tZZ[[phi.idx]] / sigma.squared + Lambda.inv[[phi.idx]])
  	b <- tZ[[phi.idx]] %*% T / sigma.squared
    alpha <- rMVN(A.tmp.chol, b)	
    rm(A.tmp.chol)
    rm(b)
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- T - Z[[phi.idx]] %*% alpha	
    sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
    rm(tmp)
    
    ##
    ## sample phi^2 for continuous time model
    ##
    
    phi.idx.star <- phi.idx + sample(( - phi.tune:phi.tune)[ - (phi.tune + 1)], 1)
    if(phi.idx.star >= 1 && phi.idx.star <= N.phi){
    	tmp.star <- T - Z[[phi.idx.star]] %*% alpha
    	mh1 <- - 1 / (2 * sigma.squared) * t(tmp.star) %*% tmp.star - 1 / 2 * t(alpha) %*% Lambda.inv[[phi.idx.star]] %*% alpha
    	rm(tmp.star)
    	tmp <- T - Z[[phi.idx]] %*% alpha
    	mh2 <- - 1 / (2 * sigma.squared) * t(tmp) %*% tmp - 1 / 2 * t(alpha) %*% Lambda.inv[[phi.idx]] %*% alpha
    	rm(tmp)
    	mh <- exp(mh1 - mh2)
    	if(mh > runif(1)){
    		phi.idx <- phi.idx.star
    		phi.accept <- phi.accept + 1 / n.mcmc
    	}
    }
       
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
    trend.save[, l] <- Z[[phi.idx]] %*% alpha
  	phi.idx.save[l] <- phi.idx
  }
	list(T.save = T.save, tau.squared.I.save = tau.squared.I.save, tau.squared.P.save = tau.squared.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, alpha.save = alpha.save, trend.save = trend.save, phi.accept = phi.accept, phi.idx.save = phi.idx.save)
}
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
  neighbor.lower <- params$neighbor.lower
  neighbor.upper <- params$neighbor.upper
  N.neighbor <- params$N.neighbor
  neighbor.model <- seq(neighbor.lower, neighbor.upper, length = N.neighbor)
  
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
  Z <- vector('list', length = N.neighbor)
  tZ <- vector('list', length = N.neighbor)
  Lambda <- vector('list', length = N.neighbor)
  Lambda.inv <- vector('list', length = N.neighbor)
  for(i in 1:N.neighbor){
    A <- make.neighborhood(neighbor.model[i], t)
    D <- diag(rowSums(A))
    Q <- D - A ## ICAR covariance matrix
    pc <- prcomp(Q)
    Z[[i]] <- pc$rotation[, 1:num.truncate]
    tZ[[i]] <- t(Z[[i]])
    Lambda[[i]] <- diag((pc$sdev^2)[1:num.truncate])
    Lambda.inv[[i]] <- diag(1 / (pc$sdev^2)[1:num.truncate])
  }
  
  I.t <- diag(rep(1, t))
  I.trunc <- diag(rep(1, num.truncate))
  ## Initial latent field
  neighbor.idx <- sample(1:N.neighbor, 1)
  alpha <- rMVN(chol(Lambda.inv[[neighbor.idx]]), rep(0, num.truncate))
  T <- Z[[neighbor.idx]] %*% alpha
  T[HI == 1] <- WI[HI == 1]  ## initialize the latent field with the observed measurements
  
  ## initialize variance values 
  tau.squared.I <- 1 / rgamma(1, alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(1, alpha.P, beta.P)
  Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
  tmp <- T - Z[[neighbor.idx]] %*% alpha	
  sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
  rm(tmp)
  
  ## initialize regression parameters
  beta.0 <- 0
  #   beta.0 <- rMVN(chol(tau.squared.P * solve(Delta.0)), rep(0, p))
  #   beta.0 <- rep(0, p)
  beta.1 <- rMVN(chol(tau.squared.P * 1 / Delta.1), 0)
  J <- rep(1, p)
  #   beta.1 <- rMVN(chol(tau.squared.P * solve(Delta.1)), rep(0, p))
  beta.1.mat <- matrix(c(1, beta.1 * J), nrow = t, ncol = p + 1, byrow = TRUE)
  beta.0.mat <- matrix(c(0, beta.0 * J), nrow = t, ncol = p + 1, byrow = TRUE)
  
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
  neighbor.accept <- 0
  neighbor.idx.save <- rep(0, n.mcmc)
  
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
    #       #     	A.tmp.vec[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (H[i, ] * c(1, beta.1))
    #       #       b[i] <- (H[i, ] * c(1, beta.1)) %*% Sigma.inv %*% (W[i, ] - H[i, ] * c(0, beta.0))
    #     }
    #     A.tmp.chol <- chol(diag(A.tmp.vec + 1 / sigma.squared))
    #     T <- rMVN(A.tmp.chol, b + Z %*% alpha / sigma.squared)
    T <- rMVN(chol((rowSums((H * beta.1.mat) %*% Sigma.inv * (H * beta.1.mat)) + 1 / sigma.squared) * I.t), (rowSums((H * beta.1.mat) %*% Sigma.inv * (W - beta.0.mat)) + Z[[neighbor.idx]] %*% alpha / sigma.squared))
    T.mat <- matrix(T, nrow = t, ncol = p, byrow = FALSE)
    
    #     rm(A.tmp.chol)
    #     rm(b)
    
    ##
    ## sample beta_0
    ##
    #     
    #     A.tmp <- 0
    #     b <- 0
    #     #     A.tmp <- matrix(0, p, p)
    #     #     b <- rep(0, p)
    #     for(i in 1:t){
    #       if(HI[i] == 1){
    #         A.tmp <- A.tmp + HP[i, ] %*% HP[i, ]
    #         #         A.tmp <- A.tmp + HP[i, ] %*% t(HP[i, ])
    #         b <- b + HP[i, ] %*% (WP[i, ] - beta.1 * HP[i, ] * T[i])
    #       }
    #     }
    #     beta.0 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.0[1])), b / tau.squared.P)
    #     #     beta.0 <- rMVN(chol(1 / tau.squared.P * A.tmp + Delta.0), b / tau.squared.P)
    #     beta.0.mat <- matrix(c(0, beta.0 * J), nrow = t, ncol = p + 1, byrow = TRUE)
    #     rm(A.tmp)    
    #     rm(b)
    
    ##
    ## sample beta_1
    ##
    #     
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
    #     beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1[1])), b / tau.squared.P)
    #     beta.1 <- rMVN(chol(1 / tau.squared.P * (A.tmp + Delta.1)), b / tau.squared.P)
    beta.1 <- rMVN(chol((sum(NP.t * T^2) + Delta.1) / tau.squared.P), sum((T.mat * HP) * (WP - HP * beta.0.mat[, - 1])) / tau.squared.P)
    beta.1.mat <- matrix(c(1, beta.1 * J), nrow = t, ncol = p + 1, byrow = TRUE)
    #     rm(A.tmp)
    #     rm(b)    
    
    ##
    ## sample tau.squared.I
    ##
    
    tau.squared.I <- 1 / rgamma(1, alpha.I + NI / 2, beta.I + sum((WI * HI - T * HI)^2) / 2)
    
    ##
    ## sample tau.squared.P
    ##
    
    tau.squared.P <- 1 / rgamma(1, alpha.P + NP / 2, beta.P + sum((WP - beta.1.mat[, - 1] * HP * T.mat - beta.0.mat[, - 1] * HP) * (WP - beta.1.mat[, - 1] * HP * T.mat - beta.0.mat[, - 1] * HP)) / 2)
    Sigma.inv <- diag(c(1 / tau.squared.I, rep(1 / tau.squared.P, p)))
    
    ##
    ## sample alpha
    ##
    
    alpha <- rMVN(chol(I.trunc / sigma.squared + Lambda.inv[[neighbor.idx]]), tZ[[neighbor.idx]] %*% T / sigma.squared)	
    talpha <- t(alpha)
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- T - Z[[neighbor.idx]] %*% alpha	
    sigma.squared <- 1 / rgamma(1, alpha.sigma + t / 2, beta.sigma + t(tmp) %*% tmp / 2)
    rm(tmp)
    
    ##
    ## sample num.neighbors for discrete time model
    ##
    
    neighbor.idx.star <- neighbor.idx + sample( (- neighbor.tune:neighbor.tune)[ - (neighbor.tune + 1)], 1)
    if(neighbor.idx.star > N.neighbor){
      neighbor.idx.star <- N.neighbor
    }
    if(neighbor.idx.star < 1){
      neighbor.idx.star <- 1
    }
    alpha.star <- rMVN(chol(I.trunc / sigma.squared + Lambda.inv[[neighbor.idx.star]]), tZ[[neighbor.idx.star]] %*% T / sigma.squared)  
    talpha.star <- t(alpha.star)
    mh1 <- - 1 / (2 * sigma.squared) * ( - 2 * t(T) %*% Z[[neighbor.idx.star]] %*% alpha.star + talpha.star %*% alpha.star) - 1 / 2 * talpha.star %*% Lambda.inv[[neighbor.idx.star]] %*% alpha.star
    mh2 <- - 1 / (2 * sigma.squared) * ( - 2 * t(T) %*% Z[[neighbor.idx]] %*% alpha + talpha %*% alpha) - 1 / 2 * talpha %*% Lambda.inv[[neighbor.idx]] %*% alpha
    mh <- exp(mh1 - mh2)
    if(mh > runif(1)){
      neighbor.idx <- neighbor.idx.star
      neighbor.accept <- neighbor.accept + 1 / n.mcmc
    }
    
    
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
    trend.save[, l]<- Z[[neighbor.idx]] %*% alpha
    neighbor.idx.save[l] <- neighbor.idx
  }
  list(T.save = T.save, tau.squared.I.save = tau.squared.I.save, tau.squared.P.save = tau.squared.P.save, beta.0.save = beta.0.save, beta.1.save = beta.1.save, sigma.squared.save = sigma.squared.save, alpha.save = alpha.save, trend.save = trend.save, neighbor.accept = neighbor.accept, neighbor.idx.save = neighbor.idx.save)
}
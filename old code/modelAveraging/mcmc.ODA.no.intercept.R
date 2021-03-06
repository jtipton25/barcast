##
## MCMC sampler for orthogonal data augmentation
##

mcmc.oda <- function(Y.list, X.o, H.list, params, epsilon = 0.001){ #Y.new, X.new, for log scoring rule
  
  ##
  ## functions and subroutines
  ##
  
  ##
  ## initialize fixed values
  ##
  
  n.mcmc <- params[[1]]
  alpha <- params[[2]]
  pi.prior <- params[[3]]  
  lambda <- params[[4]]
  
  T <- length(Y.list)
  n.o <- vector(length = T)
  for(i in 1:T){
    n.o[i] <- length(Y.list[[i]])
  }
  p <- dim(X.o)[2] - 1
  # number of augmented data points is p
  n.a <- p + 1
  I.a <- diag(n.a)  
  I.o <- diag(n.o)
  ## initialize random values
  
  gamma <- vector('list', length = T)
  for(i in 1:T){
    gamma[[i]] <- rbinom(p, 1, pi.prior)
  }
  #   lambda <- c(0, rgamma(p, alpha / 2, alpha/ 2))
#   lambda <- c(0, rep(1, p))
  Delta.gamma <- vector('list', length = T)
  for(i in 1:T){
    if(sum(gamma[[i]]) == 0){
      Delta.gamma[[i]] <- 0
    } else {
    Delta.gamma[[i]] <- diag(c(0, lambda[2:(p + 1)][which(gamma[[i]] == 1)]))  
#     Delta.gamma[[i]] <- diag(lambda[which(gamma[[i]] == 1)]) 
    }
  }
    
  ## subset based on gamma
#   X.o.gamma <- cbind(X.o[, 1], X.o[, 2:(p + 1)][, gamma == 1])
#   X.o.gamma <- vector('list', length = T)
#   for(i in 1:T){
#     X.o.gamma[[i]] <- X.o[, gamma[[i]] == 1]
#   }
#   
  ## 
  HX.o.list <- vector('list', length = T)
  tHX.o.list <- vector('list', length = T)
  for(i in 1:T){
    HX.o.list[[i]] <- X.o[H.list[[i]], ]
    tHX.o.list[[i]] <- t(HX.o.list[[i]])
  }



  ## sample Y.a
  #  Y.a <- 
  # for first sample of Y.a
  Y.a <- vector('list', length = T)
  
  ##
  ## orthogonal data augmentation
  ## 
  
#   delta <- eigen(t(X.o) %*% X.o)$values[1]
#   D <- (delta + epsilon) * I.a
#   X.a <- chol(D - (t(X.o) %*% X.o))
#   X.c <- rbind(X.o, X.a)
#   projectXontoY <- vector('list', length = T)
  
  delta <- vector('list', length = T)
  D <- vector('list', length = T)
  X.a <- vector('list', length = T)
  X.c <- vector('list', length = T)
  projectXontoY <- vector('list', length = T)

  for(i in 1:T){
    delta[[i]] <- eigen(tHX.o.list[[i]] %*% HX.o.list[[i]])$values[1]
    D[[i]] <- (delta[[i]] + epsilon) * I.a
    X.a[[i]] <- chol(D[[i]] - (tHX.o.list[[i]] %*% HX.o.list[[i]]))
    X.c[[i]] <- rbind(HX.o.list[[i]], X.a[[i]])
    projectXontoY[[i]] <- solve(t(rbind(HX.o.list[[i]], X.a[[i]])) %*% rbind(HX.o.list[[i]], X.a[[i]])) %*% t(rbind(HX.o.list[[i]], X.a[[i]]))
  }
#   X.a.gamma <- cbind(X.a[, 1], X.a[, 2:(p + 1)][, gamma == 1])
  X.a.gamma <- vector('list', length = T)
  for(i in 1:T){
    X.a.gamma[[i]] <- X.a[[i]][, gamma[[i]] == 1]
  }

  beta.tilde.gamma <- vector('list', length = T)
  for(i in 1:T){
    beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][c(TRUE, gamma[[i]] == 1), ] %*% HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] + Delta.gamma[[i]]) %*% tHX.o.list[[i]][c(TRUE, gamma[[i]] == 1), ] %*% Y.list[[i]]
  }
  
  ## initialize sigma.squared
  tmp <- vector(length = T)
  for(i in 1:T){
    tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Delta.gamma[[i]] %*% beta.tilde.gamma[[i]]
  }
  sigma.squared <- 1 / rgamma(1, sum(n.o - 1) / 2, sum(tmp) / 2)

  ## initialize variables
  O <- vector('list', length = T)
  rho <- vector('list', length = T)

  ##
  ## setup save variables
  ##
  
  gamma.save <- array(dim = c(p, T, n.mcmc))
  sigma.squared.save <- vector(length = n.mcmc)
  beta.save <- array(dim = c(p + 1, T, n.mcmc))
  rho.save <- array(dim = c(p, T, n.mcmc))
  delta.save <- delta
#   log.score.save <- vector(length = n.mcmc)  

  ##
  ## begin mcmc
  ##
  
  for(k in 1:n.mcmc){
#     if(k %% 10000 == 0){
      cat(k, ' ')
#     }
    
    ##
    ## sample Y.a
    ##
    
    for(i in 1:T){
      Y.a[[i]] <- rmvnorm(1, mean = X.a[[i]][, c(TRUE, gamma[[i]] == 1)] %*% beta.tilde.gamma[[i]], sigma = sigma.squared * (I.a + X.a[[i]][, c(TRUE, gamma[[i]] == 1)] %*% solve(tHX.o.list[[i]][c(TRUE, gamma[[i]] == 1), ] %*% HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] + Delta.gamma[[i]]) %*% t(X.a[[i]][, c(TRUE, gamma[[i]] == 1)])))
    }
    Y.c <- vector('list', length = T)
    for(i in 1:T){
      Y.c[[i]] <- c(Y.list[[i]], Y.a[[i]])
    }
    beta.hat <- vector('list', length = T)
    for(i in 1:T){
      beta.hat[[i]] <- projectXontoY[[i]] %*% Y.c[[i]]
    }
    
    ##
    ## sample sigma.squared
    ##
    
    tmp <- vector(length = T)
    for(i in 1:T){
      tmp[i] <- t(Y.list[[i]] - HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] %*% beta.tilde.gamma[[i]]) %*% (Y.list[[i]] - HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] %*% beta.tilde.gamma[[i]]) + t(beta.tilde.gamma[[i]]) %*% Delta.gamma[[i]] %*% beta.tilde.gamma[[i]]
    }
    sigma.squared <- 1 / rgamma(1, sum(n.o - 1) / 2, sum(tmp) / 2)
    
    ##
    ## sample gammma
    ##

    for(i in 1:T){
      O[[i]] <- log(pi.prior / (1 - pi.prior) * (lambda[2:(p+1)] / (delta[[i]] + lambda[2:(p+1)]))^(1 / 2)) + 1 / 2 * delta[[i]] / (delta[[i]] + lambda[2:(p+1)]) * beta.hat[[i]][2:(p + 1)]^2 / sigma.squared * delta[[i]]
      rho[[i]] <- exp(O[[i]]) / (1 + exp(O[[i]]))
      for(j in 1:p){
        if(is.na(rho[[i]][j])){
          rho[[i]][j] <- 1
        }
      }
    }
    
    for(i in 1:T){
      gamma[[i]] <- rbinom(p, 1, rho[[i]])
      if(sum(gamma[[i]]) == 0){
        Delta.gamma[[i]] <- 0
      } else {
      Delta.gamma[[i]] <- diag(lambda[c(TRUE, which(gamma[[i]] == 1))])
      }
    }
    
    ##
    ## sample beta.tilde.gamma
    ##
    
#     X.o.gamma <- cbind(X.o[, 1], X.o[, 2:(p + 1)][, gamma == 1])
#     X.a.gamma <- cbind(X.a[, 1], X.a[, 2:(p + 1)][, gamma == 1])
#     beta.tilde.gamma <- solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.o.gamma) %*% Y.o
   for(i in 1:T){
      beta.tilde.gamma[[i]] <- solve(tHX.o.list[[i]][c(TRUE, gamma[[i]] == 1), ] %*% HX.o.list[[i]][, c(TRUE, gamma[[i]] == 1)] + Delta.gamma[[i]]) %*% tHX.o.list[[i]][c(TRUE, gamma[[i]] == 1), ] %*% Y.list[[i]]
    }
    
    ##
    ## log scoring rule
    ##
    
#     log.score <- sum(dnorm(Y.new, mean = cbind(X.new[, 1], X.new[, 2:(p + 1)][, gamma == 1]) %*% beta.tilde.gamma, sd = sqrt(sigma.squared), log = TRUE))
    
    ##
    ## save samples
    ##
    
    gamma.save[, , k] <- matrix(unlist(gamma), nrow = p, ncol = T, byrow = FALSE)
    sigma.squared.save[k] <- sigma.squared
    beta.save[, , k] <- matrix(unlist(beta.hat), nrow = p + 1, ncol = T, byrow = FALSE)
    rho.save[, , k] <- matrix(unlist(rho), nrow = p, ncol = T, byrow = FALSE) 
#     log.score.save[k] <- log.score
  }
  list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, rho.save = rho.save, delta.save = delta.save)#, log.score.save = log.score.save)
  
}


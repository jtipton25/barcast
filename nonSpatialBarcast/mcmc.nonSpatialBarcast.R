##
## Model for analysis of tree ring chronologies at one spatial location ignoring different species
##

##
##
##
## Note that in the code X represents the field T in the Barcast Model
mcmc <- function(WI, W, n.mcmc, mu.0.tilde, sigma.squared.0.tilde, alpha.alpha, beta.alpha, mu.0, sigma.squared.0, alpha.sigma.squared, beta.sigma.squared, alpha.phi, beta.phi, alpha.I, beta.I, alpha.P, beta.P, mu.beta.0, sigma.squared.beta.0, mu.beta.1, sigma.squared.beta.1){
  
  ##
  ## libraries and subroutines
  ##
  
  ##
  ## initialize variables
  ##
  
  W <- cbind(WI, WP)
  t <- dim(W)[1]
  n <- dim(W)[2]
#   t.u <- c(1:451, 556) #years PDSI is unobserved
#   t.o <- 452:555 # years PDSI is observed
  HI <- !is.na(WI)
  HP <- !is.na(WP)
  H <- !is.na(W)
  X <- matrix(NA, nrow = t, ncol = n)
  X.tmp <- vector(length = t)
  NIt <- vector(length = t)
  NPt <- vector(length = t)
  for(i in 1:t){
    NIt[i] <- sum(!is.na(WI[i]))
    NPt[i] <- sum(!is.na(WP[i, ]))
  }
  MI <- sum(NIt)  
  MP <- sum(NPt)

  alpha <- runif(1, alpha.alpha, beta.alpha)
  mu <- rnorm(1, mu.0, sigma.squared.0)
  sigma.squared <- 1 / rgamma(1, alpha.sigma.squared, beta.sigma.squared)
  phi <- 1 / rgamma(1, alpha.phi, beta.phi)
  tau.squared.I <- 1 / rgamma(alpha.I, beta.I)
  tau.squared.P <- 1 / rgamma(alpha.P, beta.P)
  beta.0 <- rnorm(1, mu.beta.0, sigma.squared.beta.0)
  beta.1 <- rnorm(1, mu.beta.1, sigma.squared.beta.1)
  
  ##
  ## set up save variables
  ##
  
  X.save <- array(nrow = t, ncol = n.mcmc)
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
    
    A.chol <- chol(diag(alpha^2 / sigma.squared.epsilon * I.obs + 1))
    b <- as.numeric(by(alpha * Y[idx.na] / sigma.squared.epsilon, year.mat[idx.na], sum))
    
    X.tmp <- rMVN(A.chol, b)
    X[ - t.o] <- X.tmp[ - t.o]
    X.mat <- matrix(rep(X, n), nrow = t, ncol = n)
    
    ##
    ## sample alpha
    ##
    
    A.chol <- sqrt(sum(X^2 * I.obs / sigma.squared.epsilon + 1 / sigma.squared.alpha))
    b <- sum(X.mat[idx.na] * Y[idx.na]) / sigma.squared.epsilon + mu.alpha / sigma.squared.alpha
    alpha <- rMVN(A.chol, b)
    
    ##
    ## sample sigma.squared.epsilon
    ##
    
    #     alpha.X.tilde.mat <- matrix(rep(t(alpha) %*% t(X.tilde), n), nrow = t - 1, ncol = n)
    sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon + sum(I.obs) / 2, beta.epsilon + sum((Y[idx.na] - alpha * X.mat[idx.na])^2) / 2)
    
    ##
    ## sample mu.alpha
    ##
    
    A.chol <- sqrt(1 / sigma.squared.alpha + 1 / sigma.squared.0)
    b <- alpha / sigma.squared.alpha + mu.0 / sigma.squared.0
    mu.alpha <- rMVN(A.chol, b)
    
    ##
    ## sample sigma.squared.alpha
    ##
    
    sigma.squared.alpha <- 1 / rgamma(1, alpha.alpha + 1 / 2, beta.alpha + (alpha - mu.alpha)^2 / 2) 
    
    ##
    ## save variables
    ##
    
    X.save[, l] <- X.tmp
    alpha.save[l] <- alpha
    mu.alpha.save[l] <- mu.alpha
    sigma.squared.alpha.save[l] <- sigma.squared.alpha
    sigma.squared.epsilon.save[l] <- sigma.squared.epsilon
  }    
  
  list(X.save = X.save, alpha.save = alpha.save, mu.alpha.save = mu.alpha.save, sigma.squared.alpha.save = sigma.squared.alpha.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save)
  
}
scale.predictor <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  scale <- matrix(nrow = p, ncol = 2)
  X.tmp <- X
  for(i in 1:p){
    scale[i, ] <- c(mean(X[, i]), sqrt((n - 1) / n) * sd(X[, i]))
    X.tmp[, i] <- (X[, i] - scale[i, 1]) / scale[i, 2]
  }
  list(X = X.tmp, scale = scale)
}
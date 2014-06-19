make.output.plot <- function(out, file = NULL, method){
  if(!is.null(file)){
    jpeg(file = file, width = 6, height = 6, quality = 100, res  = 600, units = 'in')  
  }
  
  layout(matrix(1:9, 3, 3))
  matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('blue', alpha.f = 0.5))
  lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025), col = adjustcolor('red', alpha.f = 0.25))
  lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975), col = adjustcolor('red', alpha.f = 0.25))
  abline(h = 0)
  
  matplot(WI[t.o], type = 'l', col = adjustcolor('blue', alpha.f = 0.5))
  matplot(apply(out$T.save[, n.burn:n.mcmc], 1, mean)[t.o], type = 'l', add = TRUE)
  lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025)[t.o], col = adjustcolor('red', alpha.f = 0.25))
  lines(apply(out$T.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975)[t.o], col = adjustcolor('red', alpha.f = 0.25))
  abline(h = 0)
  
  matplot(out$sigma.squared.save[n.burn:n.mcmc], type = 'l', main = paste('sigma.squared = ', round(mean(out$sigma.squared.save[n.burn:n.mcmc]), digits = 4)))
  matplot(out$tau.squared.I.save[n.burn:n.mcmc], type = 'l', main = paste('tau.squared.I = ', round(mean(out$tau.squared.I.save[n.burn:n.mcmc]), digits = 4)))
  matplot(out$tau.squared.P.save[n.burn:n.mcmc], type = 'l', main = paste('tau.squared.P = ', round(mean(out$tau.squared.P.save[n.burn:n.mcmc]), digits = 4)))
  
  matplot(apply(out$trend.save[, n.burn:n.mcmc], 1, mean), type = 'l', col = adjustcolor('red', alpha.f = 0.5), ylim = c(-0.5, 0.5), main = "Trend surface")
  lines(apply(out$trend.save[, n.burn:n.mcmc], 1, quantile, prob = 0.025), type = 'l', col = adjustcolor('red', alpha.f = 0.25))
  lines(apply(out$trend.save[, n.burn:n.mcmc], 1, quantile, prob = 0.975), type = 'l', col = adjustcolor('red', alpha.f = 0.25))
  abline(h = 0)
  
  matplot(out$beta.1.save, type = 'l', main = 'Beta.1')
  
  if(method == 'continuous'){
    matplot(out$phi.idx.save, type = 'l', main = paste('accept', round(out$phi.accept / n.mcmc, digits = 4)))
  }
  
  if(!is.null(file)){
    dev.off()
  }
}
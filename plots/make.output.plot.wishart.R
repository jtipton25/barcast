make.output.plot <- function(out, resid = FALSE, file = 'filepath'){
	if(file != 'filepath'){
	  jpeg(file = file, width = 6, height = 6, quality = 100, res = 600, units = 'in')
	}
	if(resid == FALSE){
			layout(matrix(1:9, nrow = 3, ncol = 3))
			matplot(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, median), type = 'l', main = 'fitted PDSI', ylab = 'PDSI', ylim = c(-4, 4))
			steps <- 40
			for(k in (2:steps) / steps){
			  polygon(c(1:t, t:1), c(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = (k - 1 / steps) / 2), rev(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = k / 2))), col =  adjustcolor('black', alpha.f = k), border = FALSE)
        polygon(c(1:t, t:1), c(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 1 - k / 2), rev(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 1 - (k - 1 / steps) / 2))), col =  adjustcolor('black', alpha.f = k), border = FALSE)
			}
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(WI, col = adjustcolor('blue', alpha = 0.5))
			
      # plot only current period
			matplot(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, mean)[t.o], type = 'l', main = 'fitted PDSI Instrumental', ylab = 'PDSI', ylim = c(-4, 4))
			for(k in (2:steps) / steps){
			  polygon(c(1:length(t.o), length(t.o):1), c(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = (k - 1 / steps) / 2)[t.o], rev(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = k / 2)[t.o])), col =  adjustcolor('black', alpha.f = k), border = FALSE)
			  polygon(c(1:length(t.o), length(t.o):1), c(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 1 - k / 2)[t.o], rev(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 1 - (k - 1 / steps) / 2)[t.o])), col =  adjustcolor('black', alpha.f = k), border = FALSE)
			}
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 0.025)[t.o], col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[ - 1, n.burn:n.mcmc], 1, quantile, probs = 0.975)[t.o], col = adjustcolor('red', alpha = 0.25))
			lines(WI[t.o], col = adjustcolor('blue', alpha = 0.5))
			# 			plot(X, type = 'l', col = 'blue')
			MSPE <- sqrt((WI[t.o] - apply(out$X.save[ - 1, n.burn:n.mcmc], 1, mean)[t.o + 1])^2)
      plot(MSPE, type = 'l', ylab = 'RMSE for PDSI', main = paste('RMSE for PDSI ', round(mean(MSPE), 2)))
			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
      plot(out$mu.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu', ylab = 'mu')
      # 			plot(out$phi.save[n.burn:n.mcmc], type = 'l', main = 'trace for phi', ylab = 'phi')
			plot(out$tau.I.save[n.burn:n.mcmc], type = 'l', main = 'trace for tau^2_I', ylab = 'tau^2_I')
			plot(out$tau.P.save[n.burn:n.mcmc], type = 'l', main = 'trace for tau^2_P', ylab = 'tau^2_P')
			plot(out$beta.0.save[n.burn:n.mcmc], type = 'l', main = 'trace for beta_0', ylab = 'beta_0')
			plot(out$beta.1.save[n.burn:n.mcmc], type = 'l', main = 'trace for beta_1', ylab = 'beta_1')
# 			plot(out$sigma.squared.save[n.burn:n.mcmc], type = 'l', main = 'trace for sigma^2', ylab = 'sigma^2')
      if(is.null(dim(out$rho.save))){
        matplot(out$rho.save[n.burn:n.mcmc], type = 'l', main = paste('trace for rho, accept = ', round(out$rho.accept, 2)))    
      } else {
        matplot(t(out$rho.save[, n.burn:n.mcmc]), type = 'l', main = paste('trace for rho, accept = ', round(out$rho.accept, 2)))     
      }
	}
# 	} else {
# 			layout(matrix(1:6, nrow = 3, ncol = 2))
# 			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
# 			abline(h = 0, col = 'blue')
# 			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
# 			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
# 			lines(X, col = adjustcolor('blue', alpha = 0.5))
# 			# 			plot(X, type = 'l', col = 'blue')
# 			plot(sqrt((X - apply(out$X.save[, n.burn:n.mcmc], 1, mean))^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
# 			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
# 			plot(out$mu.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
# 			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
# 			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
# 	}
	if(file != 'filepath'){
		dev.off()
	}
}

##
## Plot residual diagnostics
##
# 
# make.resid.plot <- function(out, transform = 'NULL', file = 'filepath'){#model, file = 'filepath'){
# 	if(transform == 'NULL'){
# 		y.pred.out <- apply(out$y.pred.save, 1:2, mean)
# 		year.mat <- matrix(rep(1:t, n), nrow = t, ncol = n)
# 		Y.hat <- y.pred.out[idx.na]
# 		resid <- Y[idx.na] - Y.hat
# 	}
# 	if(transform == 'log'){
# 		y.pred.out <- apply(out$y.pred.save, 1:2, mean)
# 		year.mat <- matrix(rep(1:t, n), nrow = t, ncol = n)
# 		Y.hat <- y.pred.out[idx.na]
# 		resid <- log(Y[idx.na]) - Y.hat
# 	}
# 	if(file != 'filepath'){
# 		png(file)
# 	}
# 	# 	if(model == 'ar'){
# 	layout(matrix(1:4, nrow = 2, ncol = 2))
# 	plot(year.mat[idx.na], resid, main = 'residuals through time')
# 	abline(h = 0, col = 'red')
# 	plot(Y.hat, resid, main = 'residuals vs fitted')
# 	abline(h = 0, col = 'red')
# 	qqnorm(resid, main = 'Q-Q plot for residuals')
# 	qqline(resid)
# 	hist(resid, main = 'histogram of residuals')
# 	# 	}
# 	if(file != 'filepath'){
# 		dev.off()
# 	}
# }
# 

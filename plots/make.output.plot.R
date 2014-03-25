make.output.plot <- function(out, model, resid = FALSE, file = 'filepath'){
	if(file != 'filepath'){
		png(file)
	}
	if(resid == FALSE){
		if(model == 'simple'){
			layout(matrix(1:6, nrow = 3, ncol = 2))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X[t.o] ~ (1:t)[t.o], col = adjustcolor('blue', alpha = 0.5))
			plot(sqrt((X[t.o] - apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.o])^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
			plot(out$mu.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
		if(model == 'ar'){
			layout(matrix(1:9, nrow = 3, ncol = 3))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X[t.o] ~ (1:t)[t.o], col = adjustcolor('blue', alpha = 0.5))
			plot(sqrt((X[t.o] - apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.o])^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			#   polygon(c(t:1, 1:t), c(rev(apply(out$X.save[, n.burn:n.mcmc], 1, mean) - 2 * apply(out$X.save[, n.burn:n.mcmc], 1, sd)), apply(out$X.save[, n.burn:n.mcmc], 1, mean) + 2 * apply(out$X.save[, n.burn:n.mcmc], 1, sd)), col = adjustcolor('grey80', alpha.f=0.5), border = NA)  
			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
			plot(out$phi.save[n.burn:n.mcmc], type = 'l', main = paste('phi, accept = ', round(out$phi.accept, 2)), ylab = 'phi')
			plot(out$mu.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
		if(model == 'lag.two'){
			layout(matrix(1:6, nrow = 2, ncol = 3))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X, col = adjustcolor('blue', alpha = 0.5))
			plot(sqrt((X[t.o] - apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.o])^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			# 			plot(sqrt((X - apply(out$X.save[, n.burn:n.mcmc], 1, mean))^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			matplot(t(out$alpha.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for alpha', ylab = 'alpha')
			matplot(t(out$mu.alpha.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
		if(model == 'iCAR'){
			layout(matrix(1:9, nrow = 3, ncol = 3))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X, col = adjustcolor('blue', alpha = 0.5))
			# 			plot(X, type = 'l', col = 'blue')
			plot(sqrt((X[t.o] - apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.o])^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			matplot(t(out$beta.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for beta', ylab = 'beta')
			matplot(t(out$mu.beta.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for mu.beta', ylab = 'mu.beta')
			matplot(t(out$alpha.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for alpha', ylab = 'alpha')
			plot(out$sigma.squared.beta.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.beta')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
	} else {
		if(model == 'simple'){
			layout(matrix(1:6, nrow = 3, ncol = 2))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X, col = adjustcolor('blue', alpha = 0.5))
			# 			plot(X, type = 'l', col = 'blue')
			plot(sqrt((X - apply(out$X.save[, n.burn:n.mcmc], 1, mean))^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
			plot(out$mu.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
		if(model == 'ar'){
			layout(matrix(1:9, nrow = 3, ncol = 3))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X, col = adjustcolor('blue', alpha = 0.5))
			#   polygon(c(t:1, 1:t), c(rev(apply(out$X.save[, n.burn:n.mcmc], 1, mean) - 2 * apply(out$X.save[, n.burn:n.mcmc], 1, sd)), apply(out$X.save[, n.burn:n.mcmc], 1, mean) + 2 * apply(out$X.save[, n.burn:n.mcmc], 1, sd)), col = adjustcolor('grey80', alpha.f=0.5), border = NA)  
			plot(sqrt((X - apply(out$X.save[, n.burn:n.mcmc], 1, mean))^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			plot(out$alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for alpha', ylab = 'alpha')
			plot(out$phi.save[n.burn:n.mcmc], type = 'l', main = paste('phi, accept = ', round(out$phi.accept, 2)), ylab = 'phi')
			plot(out$mu.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
		if(model == 'lag.two'){
			layout(matrix(1:6, nrow = 2, ncol = 3))
			plot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI')
			abline(h = 0, col = 'blue')
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
			lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
			lines(X, col = adjustcolor('blue', alpha = 0.5))
			plot(sqrt((X - apply(out$X.save[, n.burn:n.mcmc], 1, mean))^2), type = 'l', main = 'RMSE for PDSI', ylab = 'RMSE for PDSI')
			matplot(t(out$alpha.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for alpha', ylab = 'alpha')
			matplot(t(out$mu.alpha.save[, n.burn:n.mcmc]), type = 'l', main = 'trace for mu.alpha', ylab = 'mu.alpha')
			plot(out$sigma.squared.alpha.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.alpha')
			plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l', main = 'trace plot', ylab = 'sigma.squared.epsilon')
		}
	}
	if(file != 'filepath'){
		dev.off()
	}
}

##
## Plot residual diagnostics
##

make.resid.plot <- function(out, transform = 'NULL', file = 'filepath'){#model, file = 'filepath'){
	if(transform == 'NULL'){
		y.pred.out <- apply(out$y.pred.save, 1:2, mean)
		year.mat <- matrix(rep(1:t, n), nrow = t, ncol = n)
		Y.hat <- y.pred.out[idx.na]
		resid <- Y[idx.na] - Y.hat
	}
	if(transform == 'log'){
		y.pred.out <- apply(out$y.pred.save, 1:2, mean)
		year.mat <- matrix(rep(1:t, n), nrow = t, ncol = n)
		Y.hat <- y.pred.out[idx.na]
		resid <- log(Y[idx.na]) - Y.hat
	}
	if(file != 'filepath'){
		png(file)
	}
	# 	if(model == 'ar'){
	layout(matrix(1:4, nrow = 2, ncol = 2))
	plot(year.mat[idx.na], resid, main = 'residuals through time')
	abline(h = 0, col = 'red')
	plot(Y.hat, resid, main = 'residuals vs fitted')
	abline(h = 0, col = 'red')
	qqnorm(resid, main = 'Q-Q plot for residuals')
	qqline(resid)
	hist(resid, main = 'histogram of residuals')
	# 	}
	if(file != 'filepath'){
		dev.off()
	}
}


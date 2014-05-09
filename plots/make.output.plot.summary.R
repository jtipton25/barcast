make.output.plot.summary <- function(out){
	layout(matrix(1:2, ncol = 2))
	matplot(apply(out$X.save[, n.burn:n.mcmc], 1, mean), type = 'l', main = 'fitted PDSI', ylab = 'PDSI', ylim = c(-4, 4))
	abline(h = 0, col = 'blue')
	lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025), col = adjustcolor('red', alpha = 0.25))
	lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975), col = adjustcolor('red', alpha = 0.25))
	lines(WI, col = adjustcolor('blue', alpha = 0.5))
#
matplot(apply(out$X.save[, n.burn:n.mcmc], 1, mean)[t.o], type = 'l', main = 'fitted PDSI', ylab = 'PDSI', ylim = c(-4, 4))
abline(h = 0, col = 'blue')
lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.025)[t.o], col = adjustcolor('red', alpha = 0.25))
lines(apply(out$X.save[, n.burn:n.mcmc], 1, quantile, probs = 0.975)[t.o], col = adjustcolor('red', alpha = 0.25))
lines(WI[t.o], col = adjustcolor('blue', alpha = 0.5))
}
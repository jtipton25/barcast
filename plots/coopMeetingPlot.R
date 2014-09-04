load('~/barcast/data/coopMeetingContinuousModelFit.RData')

year <- 2005 - (556:1)
# jpeg(file = '~/barcast/plots/continuousReconstruction.jpeg', width = 12, height = 6, quality = 100, res  = 600, units = 'in')
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, median), type = 'l', xaxt = 'n', xlab = 'Year', ylab = 'PDSI')
axis(1, at = 1:length(year), labels = year)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.025), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.975), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
#dev.off()

load('~/barcast/data/continuousFitOutputThu Jul 17 14:11:41 2014.RData')

year <- 2005 - (556:1)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, median), type = 'l', xaxt = 'n', xlab = 'Year', ylab = 'PDSI')
axis(1, at = 1:length(year), labels = year)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.025), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.975), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)

load('~/barcast/data/ContinuousData.RData')
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, median), type = 'l', xaxt = 'n', xlab = 'Year')
axis(1, at = 1:length(year), labels = year)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.025), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)
matplot(apply(out$T.save[, dim(out$T.save)[2] / 5 : dim(out$T.save)[2]], 1, quantile, prob = 0.975), type = 'l', col = adjustcolor('red', alpha.f = 0.25), add = TRUE)

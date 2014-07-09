##
## This model seems to be working well, maybe reduce the degrees of freedom in the Wishart matrix and run for longer time to get convergence
##

rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

library(dplR)
library(mvtnorm)
library(MASS)
suppressMessages(library(MCMCpack))
source('~/barcast/rMVN.R')
source('~/barcast/plots/make.output.plot.wishart.R')
source('~/barcast/nonSpatialBarcastWishartTakeTwo/mcmc.nonSpatialBarcast.R')

##
## load data
##

load('~/barcast/data/Hudson_5_13_2014.RData')

ls()

# jpeg(file = '~/barcast/plots/barcast_4_16_2014.jpeg', width = 6, height = 6, quality = 100, res  = 600, units = 'in')
make.output.plot(out)
# dev.off()
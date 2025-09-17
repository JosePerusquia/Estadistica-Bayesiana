####################################################################
# Normal - normal - gamma model
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(rstan)      # Version 2.32.6
library(here)       # Version 1.0.1
####################################################################

####################################################################
# Simulate some data
m=5
n=100
sigma = 10

set.seed(3141)
lambda = rgamma(1,m/2,m/2)
mu = rnorm(1,0,sqrt(1/lambda))
x = rnorm(n,mu,sqrt(sigma))
####################################################################

####################################################################
# Fit model with RStan
norm_dat =list(x=x,N=n,sigma=sqrt(sigma))

normfit = stan(file ="norm_norm_gamma.stan",
              save_dso = FALSE, iter = 1000,data=norm_dat)
print(normfit)
plot(normfit)
####################################################################
####################################################################
# Normal - Cauchy model
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
n=1
sigma = 1

set.seed(3141592)
mu = rcauchy(1)
x = rnorm(1,mu,1)
####################################################################

####################################################################
# Fit model with RStan
norm_dat =list(x=x,N=n,sigma=sqrt(sigma))

normfit = stan(file ="norm_cauchy.stan",
              save_dso = FALSE, iter = 1000,data=norm_dat)
print(normfit)
plot(normfit)
####################################################################
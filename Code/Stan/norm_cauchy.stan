//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> N;
  real x;
  real<lower=0> sigma;
}

parameters{
  real mu;
}

model {
  
  // Cauchy prior
  mu ~ cauchy(0,1);
  
  // Likelihood for observed data
  x ~ normal(mu, sigma);
}


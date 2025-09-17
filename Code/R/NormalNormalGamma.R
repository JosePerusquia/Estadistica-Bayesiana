####################################################################
# Normal - normal - gamma model
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(ggplot2)      # Version 3.5.1
library(ggthemes)     # Version 5.0.0
library(dplyr)        # Version 1.1.4
####################################################################

####################################################################
# Simulating from the marginal of mu with m=5 which results in a 
# t distribution with 5 degrees of freedom
m=5
lambda = numeric(1000)
mu_lambda = numeric(1000)

set.seed(314159)
for(i in 1:1000){
  lambda[i]=rgamma(1,m/2,m/2)
  mu_lambda[i]=rnorm(1,0,sqrt(1/lambda[i]))
}


df=data.frame(x=mu_lambda)

ggplot(data=df,aes(x=x))+
  geom_histogram(col='black',fill='skyblue2')+
  theme_minimal()+
  labs(x='',y='')
####################################################################

####################################################################
# Simulate some data
n=100
sigma = 10

set.seed(3141)
lambda = rgamma(1,m/2,m/2)
mu = rnorm(1,0,sqrt(1/lambda))
x = rnorm(n,mu,sqrt(sigma))
x_bar = mean(x)
####################################################################

####################################################################
# Posterior distribution using Gibbs
T=1000
lambda_post = numeric(T+1)
mu_post = numeric(T+1)

lambda_post[1]=1
mu_post[1]=0


for(i in 1:T){
  lambda_post[i+1] = rgamma(1,(m+1)/2,((mu_post[i]^2)+m)/2)
  sigma_1 = ((n/sigma)+lambda_post[i+1])^(-1)
  mu_1 = sigma_1*((n*x_bar)/sigma)
  mu_post[i+1] = rnorm(1,mu_1,sqrt( sigma_1))
  
}

N = c(1:(T+1))
df_post = data.frame(n=N,mu=mu_post,lambda=lambda_post)

# posterior of mu
ggplot(data=df_post,aes(x=n,y=mu))+
  geom_line(col='darkgray')+
  theme_minimal()+
  labs(x=expression(n),y=expression(mu))+
  geom_hline(yintercept = mean(mu_post),col='red')+
  geom_hline(yintercept = quantile(mu_post,.025),col='blue',
             linetype=2)+
  geom_hline(yintercept = quantile(mu_post,.975),col='blue',
             linetype=2)

mean(mu_post)  
quantile(mu_post,c(.025,.5,.975))

# posterior of lambda
ggplot(data=df_post,aes(x=n,y=lambda))+
  geom_line(col='darkgray')+
  theme_minimal()+
  labs(x=expression(n),y=expression(lambda))+
  geom_hline(yintercept = mean(lambda_post),col='red')+
  geom_hline(yintercept = quantile(lambda_post,.025),col='blue',
             linetype=2)+
  
  geom_hline(yintercept = quantile(lambda_post,.975),col='blue',
             linetype=2)

mean(lambda_post)  
quantile(lambda_post,c(.025,.5,.975))
####################################################################
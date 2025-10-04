####################################################################
# Box-Cox transform for survival times
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(ggplot2)      # Version 3.5.2
library(ggthemes)     # Version 5.1.0
library(dplyr)        # Version 1.1.4
library(MASS)         # Version 7.3-64
library(nortest)      # Version 1.0-4
####################################################################

####################################################################
# Data :survival times (from Collett, 1994) of patients in a 
# study on multiple myeloma
y=c(13,52,6,40,10,7,66,10,10,14,16,4,65,5,11,10,15,5,76,56,88,
    24,51,4,40,8,18,5,16,50,40,1,36,5,10,91,18,1,18,6,1,23,15,
    18,12,12,17,3)

df=data.frame(y)

# Skewed distribution
ggplot(data=df,aes(y))+
  geom_histogram(breaks=hist(y,plot=F)$breaks,
                 col='black',fill='skyblue4')+
  theme_minimal()+
  labs(x='',y='')
####################################################################

####################################################################
# Box-Cox transform 
w_lambda = function(y,lambda,eps = 1e-8){
  if(abs(lambda)<eps){
    w=log(y)
  }else{
    w = ((y^lambda)-1)/lambda
  }
  return(w)
}

# Box-Cox function to estimate lambda
res=boxcox(lm(y~1),plotit = T)
lambda = res$x[which.max(res$y)]

res=as.data.frame(res)
ggplot(data=res,aes(x=x,y=y))+
  geom_line()+
  labs(x=expression(lambda),y='log-likelihood')+
  theme_minimal()+
  geom_hline(yintercept=max(res$y),linetype=2)+
  geom_vline(xintercept=res$x[which(res$y==max(res$y))],linetype=2)

# Plot the transformed observations
df$w=w_lambda(y,lambda)
ggplot(data=df,aes(w))+
  geom_histogram(breaks=hist(df$w,plot=F)$breaks,
                 col='black',fill='skyblue4')+
  theme_minimal()+
  labs(x='',y='')

# Nonparametric normality tests
ad.test(df$w)
lillie.test(df$w)
shapiro.test(df$w)
####################################################################

####################################################################
# Gibbs sampling

# Samples mu from normal distribution
sample_mu = function(w,sigma2,n){
  mu = rnorm(1,mean(w),sqrt(sigma2/n))
}

# Samples sigma2 from inverse gamma distribution
sample_sigma2 = function(n,w,mu){
  rate = sum((w-mu)^2)/2
  sigma2 = rgamma(1,shape=(n-1)/2,rate=rate)
  return(1/sigma2)
}

# Log unnormalised conditional density of lambda
f_lambda=function(y,mu,sigma2,lambda){
  w=w_lambda(y,lambda)
  term1 = (lambda-1)*sum(log(y))  # Jacobian
  term2 = sum((w-mu)^2)/(2*sigma2)
  return(term1-term2)
}

# MH step with normal proposal
sample_lambda = function(lambda,mu,sigma2,y,w,tau){
  prop = rnorm(1,lambda,sqrt(tau))
  prob = f_lambda(y,mu,sigma2,prop)-f_lambda(y,mu,sigma2,lambda)
  u=runif(1)
  if(log(u)<prob){
    return(prop)
  }else{
    return(lambda)
  }
}

# Complete Gibbs algorithm
gibbs = function(y,mu_in,sigma2_in,lambda_in,tau,N){
  n=length(y)
  mu = numeric(N+1)
  sigma2 = numeric(N+1)
  lambda = numeric(N+1)
  
  mu[1]=mu_in
  sigma2[1]=sigma2_in
  lambda[1]=lambda_in
  
  for(i in 2:(N+1)){
    w=w_lambda(y,lambda[i-1])
    mu[i]=sample_mu(w,sigma2[i-1],n)
    sigma2[i]=sample_sigma2(n,w,mu[i])
    lambda[i]=sample_lambda(lambda[i-1],mu[i],sigma2[i],y,w,tau)
  }
  
  t=c(1:(N+1))
  res = data.frame(t=t,m=mu,s=sigma2,l=lambda)
  
return(res)
}

# Initial conditions
mu_in=mean(log(y))
sigma2_in=var(log(y))
lambda_in=0
tau=0.05
N=10000

# Run 10000 iterations of the algorithm
set.seed(31415)
res=gibbs(y,mu_in,sigma2_in,lambda_in,tau,N)

# Burn-in period 
B=2500

post_l = res$l[c((B+2):10001)]
post_m = res$m[c((B+2):10001)]
post_s = res$s[c((B+2):10001)]
t=c(1:7500)

df_post=data.frame(t,l=post_l,m=post_m,s=post_s)

# Plots
ggplot(df_post,aes(x=t,y=l))+
  geom_line()+
  theme_minimal()+
  labs(x='',y=expression(lambda))

ggplot(df_post,aes(x=t,y=m))+
  geom_line()+
  theme_minimal()+
  labs(x='',y=expression(mu))

ggplot(df_post,aes(x=t,y=s))+
  geom_line()+
  theme_minimal()+
  labs(x='',y=expression(sigma^2))

# Summaries
summary(df_post$l)
summary(df_post$m)
summary(df_post$s)
####################################################################
####################################################################
# Poisson change-point
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(ggplot2)      #Version 3.5.2
library(ggthemes)     #Version 5.0.0
####################################################################

####################################################################
# Data

# Sample size
n=100

# Change-point value
m=50

# Parameters of both Poisson distributions
lambda1=3
lambda2=7

# Simulate data
X1=rpois(m,lambda1)
X2=rpois(n-m,lambda2)

# Plot
x=c(1:n)
y=c(X1,X2)

df=data.frame(x=x,y=y)
ggplot(data=df,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')
####################################################################

####################################################################
# Gibbs sampling algorithm

# Function to simulate from lambda1 | lambda2, M , X
sample_l1=function(X,a1,b1,M){
  alpha_aux=a1+sum(X[c(1:M)])
  beta_aux=b1+M
  return(rgamma(1,shape=alpha_aux,rate=beta_aux))
}

# Function to simulate from lambda2 | lambda1, M, X
sample_l2=function(X,a2,b2,M){
  n=length(X)
  alpha_aux=a2+sum(X[c((M+1):n)])
  beta_aux=b2+n-M
  return(rgamma(1,shape=alpha_aux,rate=beta_aux))
}

# Function to simulate from M | lambda1, lambda2, X
sample_M=function(X,l1,l2){
  n=length(X)
  x=c(1:(n-1))
  probs=numeric(length(x))
  
  for(i in 1:length(probs)){
    probs[i]=(sum(X[c(1:i)])*log((l1/l2)))+((l2-l1)*i)
  }
  
  probs=exp(probs)
  return(sample(x,1,prob=probs))
}

# Gibbs sampler 
gibbs_sampler=function(X,N,a1,a2,b1,b2,l10,l20,M0,plot=T){
  l1=numeric(N+1)
  l2=numeric(N+1)
  Mres=numeric(N+1)
  
  l1[1]=l10
  l2[1]=l20
  Mres[1]=M0
  
  for(i in 2:(N+1)){
    l1[i]=sample_l1(X,a1,b1,Mres[i-1])
    l2[i]=sample_l2(X,a2,b2,Mres[i-1])
    Mres[i]=sample_M(X,l1[i],l2[i])
  }
  
  if(plot){
    x=c(1:N)
    y1=l1[c(2:(N+1))]
    y2=l2[c(2:(N+1))]
    y3=Mres[c(2:(N+1))]
    
    df1=data.frame(x,y1)
    p1=ggplot(data=df1,aes(x=x,y=y1))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(lambda[1]))
    plot(p1)
    
    df2=data.frame(x,y2)
    p2=ggplot(data=df2,aes(x=x,y=y2))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(lambda[2]))
    plot(p2)
    
    df3=data.frame(x,y3)
    p3=ggplot(data=df3,aes(x=x,y=y3))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(M))
    plot(p3)
  }
  L=list('L1'=l1,'L2'=l2,'M'=Mres)
  return(L)
}
####################################################################

####################################################################
# Run the algorithm with 5000 iterations and prior ditributions
# lambda1 ~ gamma(0.001,0.001), lambda2 ~ gamma(0.001,0.001) and
# M ~ Unif{1,...,100} and starting values lambda1=lambda2=1 and M=70
set.seed(314159)
res=gibbs_sampler(y,10000,.001,.001,.001,.001,1,1,70)

# Point estimates and confidence intervals
mean(res$L1)
quantile(res$L1,c(0.05,.5,.95))

mean(res$L2)
quantile(res$L2,c(0.05,.5,.95))

quantile(res$M,c(0.05,.5,.95))
####################################################################
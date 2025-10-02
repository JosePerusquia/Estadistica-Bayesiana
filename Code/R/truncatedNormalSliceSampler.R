####################################################################
# Simulating truncated normal with slice sampler
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(ggplot2)      # Version 3.5.2
library(ggthemes)     # Version 5.1.0  
####################################################################

####################################################################
# Gibbs sampling algorithm

# Function to simulate from X1
simular_X=function(Z,mu,mu_tr,sigma){
  UI=mu_tr
  US=sqrt(-2*sigma*log(Z))
  return(runif(1,UI,US))
}

# Function to simulate from Z
simular_Z=function(X,mu,sigma){
  UI=0
  US=exp(-((X-mu)^2)/(2*(sigma)))
  return(runif(1,UI,US))
}

# Gibbs sampler 
gibbs_sampler=function(N,mu,mu_tr,sigma,X0,Z0,plot=T){
  
  X_res=numeric(N+1)
  Z_res=numeric(N+1)
  X_res[1]=X0
  Z_res[1]=Z0
  
  for(i in 2:(N+1)){
    X_res[i]=simular_X(Z_res[i-1],mu,mu_tr,sigma)
    Z_res[i]=simular_Z(X_res[i],mu,sigma)
  }
  
  if(plot){
    x=c(1:N)
    y1=X_res[c(2:(N+1))]
    y2=Z_res[c(2:(N+1))]
    
    df1=data.frame(x,y1)
    p1=ggplot(data=df1,aes(x=x,y=y1))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(X))
    plot(p1)
    
    df2=data.frame(x,y2)
    p2=ggplot(data=df2,aes(x=x,y=y2))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(Z))
    plot(p2)
    
  }
  L=data.frame('X'=X_res[-1],'Z'=Z_res[-1])
  return(L)
  
}

####################################################################

####################################################################
# Density of the truncated normal  
dtruncnorm=function(x,mu,sd,a,b){
  if(x>=a &&x<=b){
    return((1/sd)*(dnorm((x-mu)/sd)/(pnorm((b-mu)/sd)-
                                       pnorm((a-mu)/sd))))
  }else{
    return(0)
  }
}

# Run the algorithm for 1000 iterations
set.seed(3141592)

mu=0
mu_tr=3
sigma=1
X0=0
Z0=exp(-((mu-mu_tr)^2)/(2*sigma))/2

res=gibbs_sampler(1000,mu,mu_tr,sigma,X0,Z0)

# Histogram of X
x=seq(mu_tr,5,by=.01)
y=numeric(length(x))

for(i in 1:length(x)){
  y[i]=dtruncnorm(x[i],mu,sigma,mu_tr,Inf)
}

df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=X,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$X,plot=F)$breaks,
                 fill='skyblue4',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()
####################################################################

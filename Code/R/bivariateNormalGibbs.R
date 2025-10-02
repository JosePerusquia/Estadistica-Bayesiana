####################################################################
# Simulating bivariate normal with Gibbs
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

# Function to sample from X1 give X2=x2
simular_X1=function(X2,mu1,mu2,sigma1,sigma2,sigma12){
  mu_aux=mu1+(sigma12/(sigma2^2))*(X2-mu2)
  sigma_aux=sqrt(sigma1^2-(sigma12/sigma2)^2)
  return(rnorm(1,mean=mu_aux,sd=sigma_aux))
}

# Function to sample from X2 give X1=x1
simular_X2=function(X1,mu1,mu2,sigma1,sigma2,sigma12){
  mu_aux=mu2+(sigma12/(sigma1^2))*(X1-mu1)
  sigma_aux=sqrt(sigma2^2-(sigma12/sigma1)^2)
  return(rnorm(1,mean=mu_aux,sd=sigma_aux))
}

# Gibbs sampler 
gibbs_sampler=function(N,mu1,mu2,sigma1,sigma2,sigma12,
                       X10,X20,plot=T){
  
  X1_res=numeric(N+1)
  X2_res=numeric(N+1)
  X1_res[1]=X10
  X2_res[1]=X20
  
  for(i in 2:(N+1)){
    X1_res[i]=simular_X1(X2_res[i-1],mu1,mu2,sigma1,
                         sigma2,sigma12)
    X2_res[i]=simular_X2(X1_res[i],mu1,mu2,sigma1,
                         sigma2,sigma12)
  }
  
  if(plot){
    x=c(1:N)
    y1=X1_res[c(2:(N+1))]
    y2=X2_res[c(2:(N+1))]

    df1=data.frame(x,y1)
    p1=ggplot(data=df1,aes(x=x,y=y1))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(X[1]))
    plot(p1)
    
    df2=data.frame(x,y2)
    p2=ggplot(data=df2,aes(x=x,y=y2))+
      geom_line()+
      theme_minimal()+
      labs(x='n',y=expression(X[2]))
    plot(p2)

  }
  L=data.frame('L1'=X1_res,'L2'=X2_res)
  return(L)
  
}

####################################################################

####################################################################
# Run the algorithm with 1000 iterations for the case
# mu1=mu2=0 , sigma1=sigma2=1 , sigma12=.3 (low correlation)
set.seed(3141592)
res=gibbs_sampler(1000,0,0,1,1,.3,0,0)

# Histogram of X1
x=density(res$L1)$x
y=density(res$L1)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L1,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$L1,plot=F)$breaks,
                 fill='skyblue4',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Histogram of X2
x=density(res$L2)$x
y=density(res$L2)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L2,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$L2,plot=F)$breaks,
                 fill='skyblue4',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Autocorrelation for X1
cor_X1=acf(res$L1,plot=F,lag.max=30)
lag = cor_X1$lag
val = cor_X1$acf
u = rep(1.96/sqrt(1000),31)
l = -u
cor_X1 = data.frame(lag,val,u,l)  

ggplot(data=cor_X1)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))

# Autocorrelation for X2
cor_X2=acf(res$L2,plot=F,lag.max=30)
lag = cor_X2$lag
val = cor_X2$acf
cor_X2 = data.frame(lag,val,u,l)  

ggplot(data=cor_X2)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))
####################################################################

####################################################################
# Run the algorithm 1000 iterations for the initial conditions
# mu1=mu2=0 , sigma1=sigma2=1 , sigma12=.99 (strong correlation)
set.seed(3141592)
res=gibbs_sampler(1000,0,0,1,1,.99,0,0)

# Histogram of X1
x=density(res$L1)$x
y=density(res$L1)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L1,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$L1,plot=F)$breaks,
                 fill='skyblue4',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Histogram of X2
x=density(res$L2)$x
y=density(res$L2)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L2,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$L2,plot=F)$breaks,
                 fill='skyblue4',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Autocorrelation for X1
cor_X1=acf(res$L1,plot=F,lag.max=30)
lag = cor_X1$lag
val = cor_X1$acf

cor_X1 = data.frame(lag,val,u,l)  

ggplot(data=cor_X1)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))

# Autocorrelation for X2
cor_X2=acf(res$L2,plot=F,lag.max=30)
lag = cor_X2$lag
val = cor_X2$acf
cor_X2 = data.frame(lag,val,u,l)  

ggplot(data=cor_X2)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))
####################################################################
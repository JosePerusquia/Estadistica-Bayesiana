#################################################################
# Importance sampling                    
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
#################################################################

#################################################################
# Libraries 
library(ggplot2)      # Version 3.5.1
library(ggthemes)     # Version 5.0.0
library(dplyr)        # Version 1.1.4
#################################################################

#################################################################
# Plots of the objective function and instrumental functions
x=seq(-5,5,by=.01)
f_x=dt(x,3)
xf_x=abs(x)*f_x
g_x=dt(x,1)
h_x=dnorm(x,0,1)

X=rep(x,4)
Y=c(xf_x,f_x,g_x,h_x)
group=c(rep('1',length(x)),rep('2',length(x)),
        rep('3',length(x)),rep('4',length(x)))

df=data.frame(X,Y,group)

ggplot(data=df,aes(x=X,y=Y,colour=group))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", 
                       labels = c(bquote("|"~X~"|"~t[3]),
                                  expression(t[3]),
                                  expression(t[1]),
                                  expression(N(0,1))))
#################################################################

#################################################################
# Monte Carlo integration with direct sampling 

# Number of samples considered and size of each sample
n=500
m=100

# Matrix to save the absolute value of the samples
sims=matrix(nrow=m,ncol=n)
set.seed(31415)
for(i in 1:m){
  sims[i,]=abs(rt(n,3))  
}

# Matrix of estimates for each n and m
est=matrix(nrow=m,ncol=n)
for(i in 1:m){
  est[i,]=cummean(sims[i,])
}

# Select the first one as the point estimate
x=c(1:n)
mu_hat=est[1,]
df_est=data.frame(x,mu_hat)

# Range of the estimator for each value of n
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

# Plot of the point estimate and the range of the estimator
X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)
df_range=data.frame(X,Y)
ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line(linewidth=.7)+
  theme_minimal()+
  labs(x='',y='')
#################################################################

#################################################################
# Importance sampling with t_1 instrumental distribution

# Matrix to save the samples
sims=matrix(nrow=m,ncol=n)

set.seed(31415)
for(i in 1:m){
  Xi=rt(n,1)
  w_x=dt(Xi,3)/dt(Xi,1)
  sims[i,]=abs(Xi)*w_x  
}

# Matrix of estimators
est=matrix(nrow=m,ncol=n)
for(i in 1:m){
  est[i,]=cummean(sims[i,])
}

# The first one as the point estimate
x=c(1:n)
mu_hat=est[1,]
df_est=data.frame(x,mu_hat)

# Range of the estimator
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

# Plot
X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)

df_range=data.frame(X,Y)
ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line(linewidth=.7)+
  theme_minimal()+
  labs(x='',y='')

# Plot of 20 weights
X=seq(-5,5,by=.01)
W=dt(X,3)/dt(X,1)
df_w=data.frame(X,W)

X=rep(0,20)
set.seed(31415)
Y=sort(rt(20,1))
W_y=dt(Y,3)/dt(Y,1)
df_y=data.frame(X,Y,W_y)

ggplot(data=df_w,aes(x=X,y=W))+
  geom_line()+
  geom_point(data=df_y,aes(x=Y,y=W_y))+
  geom_segment(data=df_y,aes(x=Y,xend=Y,y=X,yend=W_y))+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(xlim=c(-4,4))
#################################################################

#################################################################
# Importance sampling with N(0,1) instrumental distribution

# Matrix to save the samples
sims=matrix(nrow=m,ncol=n)

set.seed(31415)
for(i in 1:m){
  Xi=rnorm(n)
  w_x=dt(Xi,3)/dnorm(Xi)
  sims[i,]=abs(Xi)*w_x  
}

# Matrix of estimators
est=matrix(nrow=m,ncol=n)
for(i in 1:m){
  est[i,]=cummean(sims[i,])
}

# Select the first one as the point estimate
x=c(1:n)
mu_hat=est[1,]
df_est=data.frame(x,mu_hat)

# Range of the estimator
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

# Plot
X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)
df_range=data.frame(X,Y)

ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line(linewidth=.7)+
  theme_minimal()+
  labs(x='',y='')

# Plot for 20 points
X=seq(-5,5,by=.01)
W=dt(X,3)/dnorm(X)
df_w=data.frame(X,W)

X=rep(0,20)
set.seed(31415)
Y=sort(rnorm(20))
W_y=dt(Y,3)/dnorm(Y)
df_y=data.frame(X,Y,W_y)

ggplot(data=df_w,aes(x=X,y=W))+
  geom_line()+
  geom_point(data=df_y,aes(x=Y,y=W_y))+
  geom_segment(data=df_y,aes(x=Y,xend=Y,y=X,yend=W_y))+
  theme_minimal()+
  labs(x='',y='')+
  coord_cartesian(xlim=c(-3,3),ylim=c(0,5))
#################################################################


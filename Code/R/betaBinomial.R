####################################################################
# Beta - Binomial Model
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
####################################################################

####################################################################
# Libraries
library(ggplot2)    # Version 3.5.1
library(ggthemes)   # Version 5.0.0
library(dplyr)      # Version 1.1.4
####################################################################

####################################################################
# Parameters
n=10
m=15
alpha=1
beta=3
####################################################################

####################################################################
# Simulation
set.seed(314159)
theta = rbeta(1,alpha,beta)
x = rbinom(1,n,theta)

# Prior and posterior
theta_sop = seq(0,1,by=.01)
prior = dbeta(theta_sop,alpha,beta)
post = dbeta(theta_sop,alpha+x,beta+n-x)

sop =rep(theta_sop,2)
fx = c(prior,post)
df_plot = data.frame(sop,fx)
df_plot=df_plot%>%
  mutate(class=c(rep("Prior",101),rep("Posterior",101)))

ggplot(data=df_plot,aes(x=sop,y=fx,col=class))+
  geom_line()+
  theme_minimal()+
  labs(x="",y="")+
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10))
####################################################################

####################################################################
# Predictive distribution
theta_post = numeric(100000)
z_x = numeric(100000)

set.seed(314159)
for(i in 1:100000){
  theta_post[i]=rbeta(1,alpha+x,beta+n-x)
  z_x[i] = rbinom(1,m,theta_post[i])
}

df = data.frame(z = z_x)
ggplot(data=df,aes(x=z))+
  geom_histogram(breaks=hist(z_x,plot=F)$breaks,
                 col='black',fill='skyblue3')+
  theme_minimal()+
  labs(x="",y="")

mean(z_x)
var(z_x)
####################################################################

####################################################################
# Theoretical posterior mean and variance
post_mu = m*((alpha+x)/(alpha+beta+n));post_mu
post_var = m^2*((alpha+x)*(beta+n-x))/(((alpha+beta+n)^2)*(alpha+beta+n+1)) +
           m*(beta(alpha+x+1,beta+n-x+1)/beta(alpha+x,beta+n-x));post_var
####################################################################

####################################################################
# Uniform prior
n=10
alpha=1
beta=1

set.seed(314159)
theta = rbeta(1,alpha,beta)
x = rbinom(1,n,theta)

# Prior and posterior
theta_sop = seq(0,1,by=.01)
prior = dbeta(theta_sop,alpha,beta)
post = dbeta(theta_sop,alpha+x,beta+n-x)

sop =rep(theta_sop,2)
fx = c(prior,post)
df_plot = data.frame(sop,fx)
df_plot=df_plot%>%
  mutate(class=c(rep("Prior",101),rep("Posterior",101)))

ggplot(data=df_plot,aes(x=sop,y=fx,col=class))+
  geom_line()+
  theme_minimal()+
  labs(x="",y="")+
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10))
####################################################################



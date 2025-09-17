#################################################################
# Laplace approximation                    
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian statistics 
#################################################################

#################################################################
# Libraries
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Prior density
alpha=3
beta=1/3
x=seq(0,15,by=.01)
f=dgamma(x,shape=alpha,rate=beta)

# Simulate one value of a Poisson distribution
set.seed(3141592)
n=1
X=rpois(n,2)

# Posterior density
alpha_post=alpha+n*mean(X)
beta_post=beta+n
fpost=dgamma(x,shape=alpha_post,rate=beta_post)

# Laplace approximation
mu=(alpha+n*mean(X)-1)/(beta+n)
sigma=sqrt(alpha+n*mean(X)-1)/(n+beta)
flaplace=dnorm(x,mean=mu,sd=sigma)

# Plot
Sop=rep(x,3)
Y=c(f,fpost,flaplace)
lab=c(rep('1',length(x)),rep('2',length(x)),rep('3',length(x)))

df=data.frame(Sop,Y,lab)

ggplot(df,aes(x=Sop,y=Y,colour=lab))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", labels = c("Prior", 
                                             "Post.",'Laplace'))
#################################################################



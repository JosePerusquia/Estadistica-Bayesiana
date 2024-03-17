#################################################################
# Aproximación de Laplace
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Densidad inicial
alpha=3
beta=1/3
x=seq(0,8,by=.01)
f=dgamma(x,shape=alpha,rate=beta)

# Simulamos muestra de una Poisson
set.seed(3141592)
n=1
X=rpois(n,2)

# Densidad posterior
alpha_post=alpha+n*mean(X)
beta_post=beta+n
fpost=dgamma(x,shape=alpha_post,rate=beta_post)

# Aproximación de laplace
mu=(alpha+n*mean(X)-1)/(beta+n)
sigma=sqrt(alpha+n*mean(X)-1)/(n+beta)
flaplace=dnorm(x,mean=mu,sd=sigma)

# Data frame
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



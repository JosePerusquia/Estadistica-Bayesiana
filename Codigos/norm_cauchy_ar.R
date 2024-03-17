#################################################################
# Aceptación y rechazo Normal - Cauchy
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Soporte
x=seq(-6,6,by=.01)

# Densidad de la distribución normal y distribución Cauchy
fx=dnorm(x)
gx=dcauchy(x)

# Gráfica
X=rep(x,2)
Y=c(fx,gx)
group=c(rep('1',length(x)),rep('2',length(x)))
df=data.frame(X,Y,group)

ggplot(data=df,aes(x=X,y=Y,colour=group))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", 
                       labels = c("Normal", "Cauchy"))
#################################################################

#################################################################
# Multiplicamos por la constante M que hace que f(x) < M g(x)
M=sqrt(2*pi)*exp(-.5)
X=rep(x,2)
Y=c(fx,M*gx)
group=c(rep('1',length(x)),rep('2',length(x)))
df=data.frame(X,Y,group)

# Simulamos y ejecutamos aceptación y rechazo
set.seed(31415)
x=rcauchy(100)
u=runif(100)
type=numeric(length(x))

for(i in 1:100){
  if(u[i]<=dnorm(x[i])/(M*dcauchy(x[i]))){
    type[i]=1
  }else{
    type[i]=2
  }
}

y=u*M*dcauchy(x)
df_sim=data.frame(x,y,type)


# Gráfica
ggplot(data=df,aes(x=X,y=Y))+
  geom_line(aes(colour=group))+
  geom_point(data=df_sim,aes(x,y),shape=type,show.legend = F)+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", 
                       labels = c("Normal", "M Cauchy"))+
  xlim(c(-6,6))+
  ylim(c(0,.5))
#################################################################

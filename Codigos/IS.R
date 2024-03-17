#################################################################
# Muestreo por importancia
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
library(dplyr)
#################################################################

#################################################################
# Graficamos las densidades y la función de interés
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
                       labels = c('|x|t(3)', "t(3)",
                                  't(1)','N(0,1)'))
#################################################################

#################################################################
# Muestreo directo 

# Número de muestras
n=1500

# Número de simulaciones para obtener un rango en el estimador
m=100

# Matriz para guardar simulaciones
sims=matrix(nrow=m,ncol=n)

set.seed(31415)
for(i in 1:100){
  sims[i,]=abs(rt(n,3))  
}

# Matriz para obtener el estimador para cada n y cada m
est=matrix(nrow=m,ncol=n)

for(i in 1:100){
  est[i,]=cummean(sims[i,])
}

x=c(1:1500)
mu_hat=est[1,]
df_est=data.frame(x,mu_hat)

# Rango del estimador para cada n
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)

df_range=data.frame(X,Y)

ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')


#################################################################
# Muestreo por importancia con la distribución t(1)

# Matriz para guardar simulaciones
sims=matrix(nrow=m,ncol=n)

set.seed(31415)
for(i in 1:100){
  Xi=rt(n,1)
  w_x=dt(Xi,3)/dt(Xi,1)
  sims[i,]=abs(Xi)*w_x  
}

# Matriz de estimadores
est=matrix(nrow=m,ncol=n)

for(i in 1:100){
  est[i,]=cummean(sims[i,])
}

x=c(1:1500)
mu_hat=est[1,]
df_est=data.frame(x,mu_hat)

# Rango del estimador usando las 100 muestras
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)

df_range=data.frame(X,Y)

ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# Gráfica de los pesos para 20 puntos
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
# Muestreo por importancia usando N(0,1)

# Matriz de simulaciones
sims=matrix(nrow=m,ncol=n)

set.seed(31415)
for(i in 1:100){
  Xi=rnorm(n)
  w_x=dt(Xi,3)/dnorm(Xi)
  sims[i,]=abs(Xi)*w_x  
}

# Matriz con los estimadores
est=matrix(nrow=m,ncol=n)

for(i in 1:100){
  est[i,]=cummean(sims[i,])
}


x=c(1:1500)
mu_hat=est[1,]

df_est=data.frame(x,mu_hat)

# Rango
mu_hat_min=apply(est,2,min)
mu_hat_max=apply(est,2,max)

X=rep(x,2)
Y=c(mu_hat_min,mu_hat_max)

df_range=data.frame(X,Y)

ggplot(data=df_est,aes(x=x,y=mu_hat))+
  geom_line(data=df_range,aes(x=X,y=Y),alpha=2,col='lightblue')+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# Pesos para 20 puntos
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


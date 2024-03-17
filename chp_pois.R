#################################################################
# Punto de cambio Poisson
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Parámetros verdaderos y simulación

# Tamaño de muestra
n=100

# Punto de cambio
m=50

# Parámetros de la distribución Poisson
lambda1=3
lambda2=7

# Simulación
X1=rpois(m,lambda1)
X2=rpois(n-m,lambda2)

# Gráfica
x=c(1:n)
y=c(X1,X2)

df=data.frame(x=x,y=y)
ggplot(data=df,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')
#################################################################

#################################################################
# Algoritmo Gibbs

# Función para simular de lambda 1 
sample_l1=function(X,a1,b1,M){
  alpha_aux=a1+sum(X[c(1:M)])
  beta_aux=b1+M
  return(rgamma(1,shape=alpha_aux,rate=beta_aux))
}

# Función para simular de lambda 2
sample_l2=function(X,a2,b2,M){
  n=length(X)
  alpha_aux=a2+sum(X[c((M+1):n)])
  beta_aux=b2+n-M
  return(rgamma(1,shape=alpha_aux,rate=beta_aux))
}

# Función para simular de M
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

# Algoritmo de Gibbs 
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
#################################################################

#################################################################
# Corremos el algoritmo con 5000 iteraciones, distribuciones 
# iniciales gamma(.001,.001) para lambda 1 y lambda 2 y puntos
# iniciales de la cadena lambda1=lambda2=1 y M=70
set.seed(314159)
res=gibbs_sampler(y,5000,.001,.001,.001,.001,1,1,70)

# Estimadores puntuales
mean(res$L1)
mean(res$L2)
median(res$M)
#################################################################
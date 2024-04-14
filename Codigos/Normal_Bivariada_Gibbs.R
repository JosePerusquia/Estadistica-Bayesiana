#################################################################
# Simulación de normal bivariada con método Gibbs
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Algoritmo Gibbs

# Función para simular X1
simular_X1=function(X2,mu1,mu2,sigma1,sigma2,sigma12){
  mu_aux=mu1+(sigma12/(sigma2^2))*(X2-mu2)
  sigma_aux=sqrt(sigma1^2-(sigma12/sigma2)^2)
  return(rnorm(1,mean=mu_aux,sd=sigma_aux))
}

# Función para simular X2
simular_X2=function(X1,mu1,mu2,sigma1,sigma2,sigma12){
  mu_aux=mu2+(sigma12/(sigma1^2))*(X1-mu1)
  sigma_aux=sqrt(sigma2^2-(sigma12/sigma1)^2)
  return(rnorm(1,mean=mu_aux,sd=sigma_aux))
}

# Algoritmo de Gibbs 
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

#################################################################

#################################################################
# Corremos el algoritmo con 1000 iteraciones para el caso
# mu1=mu2=0 , sigma1=sigma2=1 , sigma12=.3 y puntos iniciales
# X10=X20=0

set.seed(3141592)
res=gibbs_sampler(1000,0,0,1,1,.3,0,0)

# Histograma de X1
x=density(res$L1)$x
y=density(res$L1)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L1,y=..density..))+
  geom_histogram(breaks=hist(res$L1,plot=F)$breaks,
                 fill='gray',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Histograma de X2
x=density(res$L2)$x
y=density(res$L2)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L2,y=..density..))+
  geom_histogram(breaks=hist(res$L2,plot=F)$breaks,
                 fill='gray',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

#################################################################

#################################################################
# Corremos el algoritmo con 1000 iteraciones para el caso
# mu1=mu2=0 , sigma1=sigma2=1 , sigma12=.3 y puntos iniciales
# X10=X20=0 (caso fuertemente correlacionado)

set.seed(3141592)
res=gibbs_sampler(10000,0,0,1,1,.99,0,0)

# Histograma de X1
x=density(res$L1)$x
y=density(res$L1)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L1,y=..density..))+
  geom_histogram(breaks=hist(res$L1,plot=F)$breaks,
                 fill='gray',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()

# Histograma de X2
x=density(res$L2)$x
y=density(res$L2)$y
df_norm=data.frame(x=x,y=y)

ggplot(data=res,aes(x=L2,y=..density..))+
  geom_histogram(breaks=hist(res$L2,plot=F)$breaks,
                 fill='gray',colour='black')+
  geom_line(data = df_norm,aes(x=x,y=y),colour='red')+
  labs(x='',y='')+
  theme_minimal()
#################################################################

#################################################################
# Correlación de la cadena de X1
cor_res=acf(res$L1,plot=F)

y=cor_res$acf[,,]
x=seq(0,length(y)-1,by=1)
y_in=rep(0,length(y))

df=data.frame(x1=x,xend=x,y1=y_in,yend=y)

ggplot(data=df)+
  geom_point(aes(x=x1,y=y))+
  geom_segment(aes(x=x1,y=y1,xend=xend,yend=yend))+
  geom_hline(yintercept = qnorm((1 + .95)/2)/sqrt(1001),
             linetype=2,col='blue')+
  geom_hline(yintercept = -qnorm((1 + .95)/2)/sqrt(1001),
             linetype=2,col='blue')+
  geom_hline(yintercept = 0,linewidth=.1)+
  labs(x='',y='')+
  theme_minimal()
#################################################################
#################################################################
# Metropolis - Hastings Simular Distribución Normal
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# Función para generar la cadena
MH_norm=function(N,sigma,x_in=0,plot=T){
  
  X=numeric(N+1)
  X[1]=x_in
  
  Y=numeric(N)
  status=numeric(N)
  
  for(i in 2:(N+1)){
    eps=rnorm(1,0,sd=sigma)
    Y[i-1]=X[i-1]+eps
    
    alpha=min(dnorm(Y[i-1])/dnorm(X[i-1]),1)
    u=runif(1)
    
    if(u<alpha){
      X[i]=Y[i-1]
      status[i-1]=1
    }else{
      X[i]=X[i-1]
    }
  }
  
  n=seq(1,N,by=1)
  X=X[-1]
  df_X=data.frame(n,X)
  
  n_Y=n[which(status==0)]
  Y_aux=Y[which(status==0)]
  df_Y=data.frame(n_Y,Y_aux)
  
  p=ggplot(data=df_X,aes(x=n,y=X))+
    geom_line()+
    theme_minimal()+
    geom_point(data=df_Y,aes(x=n_Y,y=Y_aux),col='blue',fill='white',shape = 21)+
    labs(x=expression(n),y=expression(X[n]),title='')+
    coord_cartesian(ylim=c(-6,6))
  
  if(plot){
    print(p)
  }  
  
  cor=acf(X,plot=F)$acf[2,,]
  
  L=list(X=X,Y=Y,alpha=mean(status),cor=cor)
  return(L)
}
#################################################################

#################################################################
# Cadenas para diferentes valores de sigma
sigmas=c(.1,1,2.38,10)

set.seed(314159)
res=MH_norm(1000,sigmas[1])
res=MH_norm(1000,sigmas[2])
res=MH_norm(1000,sigmas[3])
res=MH_norm(1000,sigmas[4])

# A partir de 100 muestras se obtiene la probabilidad de aceptar
# y la correlación entre X(t) y X(t-1)
alphas=numeric(4)
cors=numeric(4)

set.seed(314159)
for(i in 1:4){
  for(j in 1:100){
    res=MH_norm(1000,sigmas[i],plot=F)
    alphas[i]=alphas[i]+res$alpha
    cors[i]=cors[i]+res$cor
  }
  alphas[i]=alphas[i]/100
  cors[i]=cors[i]/100
}

alphas
cors
#################################################################


#################################################################
# Caminata aleatoria adaptable
mh_arw_norm=function(N,x_in,sigma_in,tau,e1,A1,plot_x=T,
                     plot_sigma=T){
  X=numeric(N+1)
  X[1]=x_in
  
  sigma=numeric(N+1)
  sigma[1]=sigma_in

  Y=numeric(N)
  status=numeric(N)
  
  for(i in 2:(N+1)){
    
    Y[i-1]=rnorm(1,sd=sigma[i-1])
    
    alpha=min(dnorm(Y[i-1])/dnorm(X[i-1]),1)
    u=runif(1)
    
    if(u<alpha){
      X[i]=Y[i-1]
      status[i-1]=1
    }else{
      X[i]=X[i-1]
    }
    
    sigma[i]=sigma[i-1]+(sigma_in/(i-1))*(alpha-tau)
  }
  
  if(plot_sigma){
    n=c(1:(N+1))
    df_sigma=data.frame(n,sigma)
    p=ggplot(data=df_sigma,aes(x=n,y=sigma))+
      geom_line()+
      theme_minimal()+
      labs(x=expression(n),y=expression(sigma[n]))+
      geom_hline(yintercept=2.38,col='red')
    
    plot(p)
  }
  
  if(plot_x){
    n=c(1:N)
    X_aux=X[-1]
    df_X=data.frame(n,X_aux)
    
    n_Y=n[which(status==0)]
    Y_aux=Y[which(status==0)]
    df_Y=data.frame(n_Y,Y_aux)
    
    q=ggplot(data=df_X,aes(x=n,y=X_aux))+
      geom_line()+
      theme_minimal()+
      geom_point(data=df_Y,aes(x=n_Y,y=Y_aux),col='blue',fill='white',shape = 21)+
      labs(x=expression(n),y=expression(X[n]),title='')+
      coord_cartesian(ylim=c(-6,6))
    plot(q)
  }
  
  cor=acf(X[-1],plot=F)$acf[2,,]
  
  L=list(X=X[-1],sigma=sigma,Y=Y,alpha=mean(status),cor=cor)
  return(L)
}

#################################################################

#################################################################
# Comparación del algoritmo MH con caminata aleatoria con sigma
# 2.38 y la caminata aleatoria adaptable iniciando en 10
# con tau=.4434, epsilon1=.0001, A1=1000

set.seed(314159)
res=MH_norm(1000,2.38)
res1=mh_arw_norm(1000,0,10,.4434,.0001,1000)


# A partir de 100 muestras se obtiene la probabilidad de aceptar
# y la correlación entre X(t) y X(t-1)
alpha=0
cors1=0

set.seed(314159)
for(j in 1:100){
    res1=mh_arw_norm(1000,0,10,.4434,.0001,1000,plot_sigma = F,
                     plot_x=F)
    alpha=alpha+res$alpha
    cors1=cors1+res$cor
}
alpha=alpha/100;alpha
cors1=cors1/100;cors1
#################################################################

#################################################################
# Muestreo adaptable de aceptación y rechazo distribución gamma
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
#################################################################

#################################################################
# La log densidad
h=function(x,a,b){
  return(((a-1)*log(x))-(x*b))
}
#################################################################

#################################################################
# La derivada de la log densidad
dh=function(x,a,b){
  return(((a-1)/x)-b)
}
#################################################################

#################################################################
# Función que actualiza los puntos tangentes
Z=function(Tk,a,b){
  z=0
  
  for (i in 1:(length(Tk)-1)){
    aux=(h(Tk[i+1],a,b)-h(Tk[i],a,b)-Tk[i+1]*dh(Tk[i+1],a,b)+
           Tk[i]*dh(Tk[i],a,b))/(dh(Tk[i],a,b)-dh(Tk[i+1],a,b))
    z=c(z,aux)
  }
  
  z=c(z,Inf)
}
#################################################################

#################################################################
# Cubierta superior
Uk=function(Tk,a,b){
  uk=c()
  for(i in 1:length(Tk)){
    aux1=h(Tk[i],a,b)
    aux2=Tk[i]
    aux3=dh(Tk[i],a,b)
    aux=substitute(expression(a+(x-b)*c),
                   list(a=aux1,b=aux2,c=aux3))
    uk=c(uk,aux)
  }
  return(uk)
}
#################################################################

#################################################################
# Función que grafica la log-densidad y la cubierta superior
ggplot_envelope=function(a,b,z,uk,Tk,Xup){
  x=seq(0,Xup,by=.01)
  f_x=dgamma(x,shape=a,rate=b)
  g_x=(x^(a-1))*exp(-x*b)
  h_x=log(g_x)
  df=data.frame(x,h_x)
  
  h_Tk=h(Tk,a,b)
  df_Tk=data.frame(Tk,h_Tk)
  
  p=ggplot(data=df,aes(x=x,y=h_x))+
    geom_line()+
    geom_point(data=df_Tk,aes(x=Tk,y=h_Tk),shape=1)+
    geom_vline(xintercept = z[c(2:(length(z)-1))],linetype=2)+
    theme_minimal()+
    labs(x='',y='')
    
  
  for(i in 1:(length(uk)-1)){
    df_aux=data.frame(x1=z[i],
                      x2=z[i+1],
                      y1=eval(uk[[i]][[2]],list(x=z[i])),
                      y2=eval(uk[[i]][[2]],list(x=z[i+1])))
    p=p+geom_segment(data=df_aux,aes(x=x1,
                                     y=y1,
                                     xend=x2,
                                     yend=y2),colour='red')
  }
  
  
  df_aux=data.frame(x1=z[length(z)-1],
                    x2=Xup,
                    y1=eval(uk[[length(z)-1]][[2]],
                            list(x=z[length(z)-1])),
                    y2=eval(uk[[length(z)-1]][[2]],
                            list(x=Xup)))
  
  p=p+geom_segment(data=df_aux,aes(x=x1,
                                   y=y1,
                                   xend=x2,
                                   yend=y2),colour='red')
  
  plot(p)
  
}
#################################################################

#################################################################
# Función para graficar densidad y función de rechazo 
ggplot_rejection=function(a,b,z,sk,C,Xup){
  x=seq(0,Xup,by=.01)
  g_x=(x^(a-1))*exp(-x*b)
  
  df=data.frame(x,g_x)
  
  f_z=numeric(length(z))
  
  for(i in 1:(length(sk)-1)){
    f_z[i]=C*eval(sk[[i]][[2]],list(x=z[i]))
  }
  
  f_z[length(sk)]=C*eval(sk[[length(z)-1]][[2]],
                         list(x=z[length(sk)]))
  df_z=data.frame(z,f_z)
  
  p=ggplot(data=df,aes(x=x,y=g_x))+
    geom_line()+
    geom_point(data=df_z,aes(x=z,y=f_z))+
    theme_minimal()+
    labs(x='',y='')
  
  for(i in 1:(length(sk)-1)){
    aux=seq(z[i],z[i+1],by=.01)
    aux2=C*eval(sk[[i]][[2]],list(x=aux))
    df_aux=data.frame(aux,aux2)
    p=p+geom_line(data=df_aux,aes(x=aux,
                                  y=aux2),colour='red')
  }

  aux=seq(z[length(z)-1],Xup,by=.01)
  aux2=C*eval(sk[[length(z)-1]][[2]],list(x=aux))
  df_aux=data.frame(aux,aux2)
  p=p+geom_line(data=df_aux,aes(x=aux,
                                y=aux2),colour='red')
  
  plot(p)
}
#################################################################

#################################################################
# Función que obtiene la cubierta para simular

Sk=function(uk,z){
  sk=c()
  C=0
  for(i in 1:length(uk)){
    f <- function(x){ exp(eval( uk[[i]][[2]]) ) }
    c=integrate(f,z[i],z[i+1])
    C=C+c[1]$value
  }
  for(i in 1:length(uk)){
    aux=substitute(expression(exp(u)/k),list(u=uk[[i]][[2]],k=C))
    sk=c(sk,aux)
  }
  l=list(sk=sk,C=C)
  return(l)
}
#################################################################

#################################################################
# Función que obtiene la cubierta inferior
Lk=function(Tk,a,b){
  lk=expression(-Inf)
  for(i in 1:(length(Tk)-1)){
    aux1=Tk[i+1]
    aux2=h(Tk[i],a,b)
    aux3=Tk[i]
    aux4=h(Tk[i+1],a,b)
    aux=substitute(expression(((a-x)*b+(x-c)*d)/a-c),
                   list(a=aux1,b=aux2,c=aux3,d=aux4))
    lk=c(lk,aux)
  }
  lk=c(lk,expression(-Inf))
  return(lk)
}
#################################################################

#################################################################
# Función para simular de la cubierta superior
sample_sk=function(sk,z,Tk,C,a,b){
  CDF=c(0)
  for(i in 1:(length(z)-1)){
    f <- function(x){ (eval(sk[[i]][[2]]) ) }
    c=integrate(f,z[i],z[i+1])
    CDF=c(CDF,c[1]$value)
  }
  CDF=cumsum(CDF)
  U=runif(1,0,1)
  n=findInterval(U,CDF)
  a1=dh(Tk[n],a,b)
  X=log((((U-CDF[n])*C*a1)/(exp(h(Tk[n],a,b)-Tk[n]*a1)))+
          exp(z[n]*a1))/a1
  return(X)
  
}
#################################################################

#################################################################
# ARS para distribución gamma de parámetros a=>1, b>0 
ARS_gamma=function(N,Tk,a,b,Xup,plot=T){
  z=Z(Tk,a,b)
  uk=Uk(Tk,a,b)
  sk=Sk(uk,z)
  lk=Lk(Tk,a,b)
  if(plot==T){
    ggplot_envelope(a,b,z,uk,Tk,Xup)
    ggplot_rejection(a,b,z,sk$sk,sk$C,Xup)
  }
  res=numeric(N)
  count=0
  for(i in 1:N){
    flag=F
    while(flag==F){
      S=sample_sk(sk$sk,z,Tk,sk$C,a,b)
      count=count+1
      W=runif(1,0,1)
      a1=findInterval(S,Tk)+1
      a2=findInterval(S,z)
      
      if(a1==1 ||a1==length(lk)){
        aux=0
      }else{
        aux=exp(eval(lk[[a1]][[2]],
                     list(x=S))-eval(uk[[a2]][[2]],list(x=S)))
      }
      
      if(W<=aux){
        res[i]=S
        flag=T
        #print('acceptance')
      }else{
        if(W<= exp(h(S,a,b)-eval(uk[[a2]][[2]],list(x=S)))){
          res[i]=S
          flag=T
          #Tk=sort(c(Tk,S))
          #z=Z(Tk,h,dh)
          #uk=Uk(Tk,h,dh)
          #sk=Sk(uk,z)
          #lk=Lk(Tk,h,dh)
          #print('acceptance')
        }else{
          Tk=sort(c(Tk,S))
          z=Z(Tk,a,b)
          uk=Uk(Tk,a,b)
          sk=Sk(uk,z)
          lk=Lk(Tk,a,b)
          print('Fail')
          #ggplot_envelope(a,b,z,uk,Tk,Xup)
          #ggplot_rejection(a,b,z,sk$sk,sk$C,Xup)
        }
      }
    }
  }
ggplot_envelope(a,b,z,uk,Tk,Xup)
ggplot_rejection(a,b,z,sk$sk,sk$C,Xup)
return(data.frame(res))
  
}
#################################################################

#################################################################
# Ejemplo para gamma (3,1)

# Cubierta inicial con 3 puntos
Tk=c(.1,1.5,12)

set.seed(314159)
res=ARS_gamma(200,Tk,3,1,20,plot=T)

# Histograma
x=seq(0,10,by=.01)
fx=dgamma(x,3,1)

df_t=data.frame(x,fx)

ggplot(data=res,aes(x=res,y=after_stat(density)))+
  geom_histogram(breaks=hist(res$res,plot=F)$breaks,
                 color="black", fill="lightblue")+
  geom_line(data=df_t,aes(x=x,y=fx),colour='red')+
  theme_minimal()+
  labs(x='',y='')
#################################################################



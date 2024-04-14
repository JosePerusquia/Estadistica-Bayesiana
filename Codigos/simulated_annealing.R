#################################################################
# Simulated annealing
# Autor: José A. Perusquía Cortés
#################################################################

#################################################################
# Librerías
library(ggplot2)
library(ggthemes)
library(plot3D)
#################################################################

#################################################################
# Función a minimizar
h=function(x,y){
  return((4-2.1*(x^2)+(x^4)/3)*x^2+x*y+4*((y^2)-1)*(y^2))
}

# Graficamos contornos 
x=sort(rep(seq(-3,3,length.out=200),200))
y=rep(seq(-2,2,length.out=200),200)

z=numeric(200*200)

for(i in 1:40000){
  z[i]=h(x[i],y[i])
}

df=data.frame(x,y,z)
ggplot(data=df,aes(x=x,y=y,z=z))+
  geom_contour(aes(colour = after_stat(level)),show.legend = F)+
  theme_minimal()+
  labs(x='',y='')

# Función que implementa el algoritmo
simulated_annealing=function(N,x_in,y_in,beta0,alpha,temp,plot=F){
  x=numeric(N+1)
  y=numeric(N+1)
  x[1]=x_in
  y[1]=y_in
  
  for(i in 2:(N+1)){
    if(temp=="geom"){
      beta_new=(alpha^(i-1)*beta0)
    }else{
      beta_new=log(1+i)/beta0
    }
    
    # Simulamos valores propuestos por una uniforme en el
    # soporte de la función
    aux1=runif(1,-3,3)
    aux2=runif(1,-2,2)
    
    # Probabilidad de aceptación
    prob=min(1,exp(-beta_new*(h(aux1,aux2)-h(x[i-1],y[i-1]))))
    u=runif(1)
    if(u<prob){
      x[i]=aux1
      y[i]=aux2
    }else{
      x[i]=x[i-1]
      y[i]=y[i-1]
    }
  }
  
  # Evaluamos la función 
  hx=h(x,y)
  
  if(plot==T){
    
    n=c(1:(N+1))
    
    df_x=data.frame(n,x)
    df_y=data.frame(n,y)
    df_h=data.frame(n,hx)
    
    p=ggplot(data=df_x,aes(x=n,y=x))+
      geom_line()+
      theme_minimal()+
      labs(x=expression(n),y=expression(x[n]),title='')
    
    plot(p)
    
    q=ggplot(data=df_y,aes(x=n,y=y))+
      geom_line()+
      theme_minimal()+
      labs(x=expression(n),y=expression(y[n]),title='')
    
    plot(q)
    
    r=ggplot(data=df_h,aes(x=n,y=hx))+
      geom_line()+
      theme_minimal()+
      labs(x=expression(n),y=expression(h[n]),title='')
    
    plot(r)

  }
  
  
  l=list(x=x,y=y,hx=hx)
  return(l)
}
#################################################################

#################################################################
# Ejemplo para una sola corrida del algoritmo
N=1000
set.seed(314159)
res=simulated_annealing(N,0,0,.5,1.1,temp="geom",plot=T)

# El valor mínimo 
res$hx[N+1]

# En los puntos
res$x[N+1]
res$y[N+1]


# Corremos el algoritmo 5000 veces 
X=numeric(5000)
Y=numeric(5000)
H=numeric(5000)

set.seed(31415)
for(i in 1:5000){
  res=simulated_annealing(N,0,0,.5,1.1,temp="geom")
  X[i]=res$x[N+1]
  Y[i]=res$y[N+1]
  H[i]=res$hx[N+1]
}

# El mínimo de la función
min(H)

# Puntos donde se tiene el mínimo
X[which(H==min(H))]
Y[which(H==min(H))]
#################################################################

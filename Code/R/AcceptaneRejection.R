#################################################################
# Normal acceptance-rejection algorithm with Cauchy proposal                
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian Statistics 
#################################################################

#################################################################
# Libraries
library(ggplot2)        # Version 3.5.1
library(ggthemes)       # Version 5.0.0
library(nortest)        # Version 1.0-4
library(dplyr)          # Version 1.1.4
#################################################################

#################################################################
# Plots of normal and Cauchy densities

sop=seq(-6,6,by=.01)
fx=dnorm(sop)
gx=dcauchy(sop)

# Dataframe
x=rep(sop,2)
y=c(fx,gx)
group=c(rep('1',length(sop)),rep('2',length(sop)))
df=data.frame(x,y,group)

# Plot
ggplot(data=df,aes(x=x,y=y,colour=group))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", 
                       labels = c("Normal", "Cauchy"))

# Consider M such that f(x) < M g(x) for all x in sop(f)
M=sqrt(2*pi)*exp(-.5)
y=c(fx,M*gx)
df=data.frame(x,y,group)

ggplot(data=df,aes(x=x,y=y,colour=group))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')+
  scale_color_discrete(name = "", 
                       labels = c("Normal", "Cauchy"))
#################################################################

#################################################################
# Acceptance-rejection algorithm
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

# Plot
p1=ggplot(data=df,aes(x=x,y=y))+
    geom_line(aes(colour=group))+
    geom_point(data=df_sim,aes(x,y),shape=type,
               show.legend = F)+
    theme_minimal()+
    labs(x='',y='')+
    scale_color_discrete(name = "", 
                       labels = c("Normal", "M Cauchy"))
p1

p1+coord_cartesian(x=c(-6,6))

# Number of accepted values
ar = df_sim%>%
  select(type)%>%
  group_by(type)%>%
  count()%>%
  mutate(acc=case_when(type==1~'Accepted',
                       type==2~'Rejected'))

ggplot(data=ar,aes(x=acc,y=n,fill=acc))+
  geom_col(col='black',show.legend = F)+
  theme_minimal()+
  labs(x='',y='')

# qqplot 
val_norm = df_sim%>%
  filter(type==1)

ggplot(data=val_norm,aes(sample=x))+
  geom_qq()+
  geom_qq_line(col='red')+
  theme_minimal()+
  labs(x='',y='')

# Normality tests
ad.test(val_norm$x)
shapiro.test(val_norm$x)
lillie.test(val_norm$x)
#################################################################

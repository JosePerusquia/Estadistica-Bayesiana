#################################################################
# Monte Carlo approximation Normal - Cauchy model               
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Bayesian Statistics 
#################################################################

#################################################################
# Libraries
library(ggplot2)
library(ggthemes)
library(dplyr)
#################################################################

#################################################################
# Initial conditions
set.seed(3141592)
theta0 = rcauchy(1)
x = rnorm(1,theta0)
#################################################################

#################################################################
# Monte Carlo for the posterior expected value for varying n
n = c(1:1000)
post_mean = numeric(length(n))

set.seed(31415)
for(i in 1:1000){
  theta = rnorm(n[i],x,1)
  numerator = sum(theta/(1+theta^2))
  denominator = sum(1/(1+theta^2))
  post_mean[i] = numerator/denominator
}

df = data.frame(n,post_mean)

ggplot(data=df,aes(n,post_mean))+
  geom_point(size=.5)+
  geom_line(linewidth=.1)+
  labs(x=expression(n),y=expression(widehat(E)(theta ~ "|" ~ X)))+
  theme_minimal()
#################################################################

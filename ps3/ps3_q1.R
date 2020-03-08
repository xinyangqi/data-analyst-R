# Stats 506, Fall 2019
# This script contains solutions for
# problem set 3, question 1.
#
# This script computes four kinds of CI 
# for the ratio of expectations of vector x and y
# The four CIs are jackknife CI, percentile CI, 
# basic bootstrap CI and normal approximation CI
#
# Author: Xinyang Qi
# Updated: Nov 7, 2019

library(datasets)
library(tidyverse)

# a
jack_CI = function(x,y,alpha=0.05){
  # Compute CI for the ratio of expectations of vector x and y,
  # based on jackknife estimate for the standard error.
  # Agr:
  #     x: vector
  #     y: vector
  #     alpha: the significance of CI.
  # Output: jackknife CI for the ratio of expectations
  #         of vector x and y
  n1=length(x)
  n2=length(y)
  m1=mean(x)
  m2=mean(y)
  m=m1/m2
  x=rep(x,n1)
  y=rep(y,n2)
  dim(x)=c(n1,n1)
  dim(y)=c(n2,n2)
  x=x-diag(diag(x))
  y=y-diag(diag(y))
  theta_i=c((colMeans(x)*n1/(n1-1))/m2,
            m1/(colMeans(y)*n2/(n2-1)))
  theta_bar=sum(theta_i)/(n1+n2)
  sigma_jack=sqrt((n1+n2-1)/(n1+n2))*
    sqrt(sum((theta_i-theta_bar)^2))
  lcb = round(m - qnorm(1-alpha/2) * sigma_jack,5)
  ucb = round(m + qnorm(1-alpha/2) * sigma_jack,5)
  CI = paste0('(',paste(lcb,ucb,sep = ','),')')
  return(CI)
}

# b
compute_CI2 = function(x,y,B=1e3,alpha=0.05){
  # Compute CI for the ratio of expectations of vector x and y,
  # based on percentile method, basic bootstrap and
  # the normal approximation with bootstrap standard error.
  # Agr:
  #     x: vector
  #     y: vector
  #     B: the multiple of bootstrap
  #     alpha: the significance of CI.
  # Output: Three kinds of CI for the ratio of expectations
  #         of vector x and y
  n1=length(x)
  n2=length(y)
  m=mean(x)/mean(y)
  # bootstrap
  boot_x=sample(x,n1*B,replace = TRUE)
  boot_y=sample(y,n2*B,replace = TRUE)
  dim(boot_x)=c(n1,B)
  dim(boot_y)=c(n2,B)
  boot_theta=colMeans(boot_x)/colMeans(boot_y)
  # compute CIs
  per_lwr=round(quantile(boot_theta,0.025),5)
  per_upr=round(quantile(boot_theta,0.975),5)
  se=sqrt(mean((boot_theta-mean(boot_theta))^2)*B/(B-1))
  basic_lwr=round(2*m-per_upr,5)
  basic_upr=round(2*m-per_lwr,5)
  normal_lwr=round(m-qnorm(1-alpha/2)*se,5)
  normal_upr=round(m+qnorm(1-alpha/2)*se,5)
  CIs=c(paste0('(',paste(per_lwr,per_upr,collapse = ','),')'),
        paste0('(',paste(basic_lwr,basic_upr,collapse = ','),')'),
        paste0('(',paste(normal_lwr,normal_upr,collapse = ','),')'))
  names(CIs)=c('per_CI','basic_CI','normal_CI')
  return(CIs)
}

# c
# Based on the ToothGrowth data in R¡¯s datasets package,
# Using each of the methods above, 
# compute a point estimate and 95% confidence intervals 
# for the ratio comparing mean odontoblast length 
# for the ¡°OJ¡± supplement type to the ¡°VC¡± supplement type
# for each level of the ¡°dose¡± variable 
# and present the result in a nicely formatted table.
data0=force(ToothGrowth)
dose=c(0.5,1.0,2.0)
x1=data0$len[data0$supp=='OJ'&data0$dose==0.5]
x2=data0$len[data0$supp=='OJ'&data0$dose==1.0]
x3=data0$len[data0$supp=='OJ'&data0$dose==2.0]
y1=data0$len[data0$supp=='VC'&data0$dose==0.5]
y2=data0$len[data0$supp=='VC'&data0$dose==1.0]
y3=data0$len[data0$supp=='VC'&data0$dose==2.0]
mean_ratios=c(mean(x1)/mean(y1),mean(x2)/mean(y2),
              mean(x3)/mean(y3))
jack_CIs=c(jack_CI(x1,y1),jack_CI(x2,y2),jack_CI(x3,y3))
per_CIs=c(compute_CI2(x1,y1)[1],compute_CI2(x2,y2)[1],
          compute_CI2(x3,y3)[1])
basic_CIs=c(compute_CI2(x1,y1)[2],compute_CI2(x2,y2)[2],
            compute_CI2(x3,y3)[2])
normal_CIs=c(compute_CI2(x1,y1)[3],compute_CI2(x2,y2)[3],
             compute_CI2(x3,y3)[3])
result=data.frame(dose,mean_ratios,jack_CIs,per_CIs,
                  basic_CIs,normal_CIs)

knitr::kable(result,digits = 3)

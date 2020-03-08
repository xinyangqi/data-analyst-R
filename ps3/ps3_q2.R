# Stats 506, Fall 2019
# This script contains solutions for
# problem set 3, question 2.
#
# This script carries out a Monte Carlo investigation 
# to computes four kinds of CI 
# for the ratio of expectations of matrices x and y
# The four CIs are jackknife CI, percentile CI, 
# basic bootstrap CI and normal approximation CI
#
# Author: Xinyang Qi
# Updated: Nov 7, 2019

library(tidyverse)

jack_repli_CI=function(x,y,alpha=0.05){
  # Compute CI for the ratio of expectations 
  # for each replicate of matrices x and y,
  # based on jackknife estimate for the standard error.
  # Agr:
  #     x: matrix
  #     y: matrix
  #     alpha: the significance of CI.
  # Output: jackknife CI for the ratio of expectations
  #         for each replicate of matrices x and y.
  stopifnot((nrow(x)==nrow(y))|ncol(x)==ncol(y))
  if (nrow(x)!=nrow(y)){
    x=t(x)
    y=t(y)
  }
  m1=rowMeans(x)
  m2=rowMeans(y)
  m=m1/m2
  x_boot=x[rep(1:nrow(x),ncol(x)),1:ncol(x)]
  y_boot=y[rep(1:nrow(y),ncol(y)),1:ncol(y)]
  x1=x
  y1=y
  dim(x1)=c(nrow(x1)*ncol(x1),1)
  dim(y1)=c(nrow(y1)*ncol(y1),1)
  x1=as.numeric(x1)
  y1=as.numeric(y1)
  theta_i=c((rowSums(x_boot)-x1)/(m2*(ncol(x)-1)),
            m1/((rowSums(y_boot)-y1)/(ncol(y)-1)))
  dim(theta_i)=c(nrow(x),(ncol(x)+ncol(y)))
  sigma_jack=sqrt((ncol(x)+ncol(y)-1)/(ncol(x)+ncol(y)))*
    sqrt(rowMeans((theta_i-m)^2)*(ncol(x)+ncol(y)))
  jack_lwr=round(m-qnorm(1-alpha/2)*sigma_jack,5)
  jack_upr=round(m+qnorm(1-alpha/2)*sigma_jack,5)
  CI=cbind(jack_lwr,jack_upr)
  return(CI)
}

# b
boot_CI=function(x,y,B=1e3,alpha=0.05){
  # Compute CI for the ratio of expectations 
  # for each replicate of matrices x and y,
  # based on percentile method, basic bootstrap and
  # the normal approximation with bootstrap standard error.
  # Agr:
  #     x: matrix
  #     y: matrix
  #     B: the multiple of bootstrap
  #     alpha: the significance of CI.
  # Output: Three kinds of CI for the ratio of expectations
  #         for each replicate of matrices x and y.
  stopifnot((nrow(x)==nrow(y))|ncol(x)==ncol(y))
  if (ncol(x)!=ncol(y)){
    x=t(x)
    y=t(y)
  }
  theta_hat=colMeans(x)/colMeans(y)
  # bootstrap
  x_boot=x[sample(1:nrow(x),nrow(x)*B,replace = TRUE),1:ncol(x)]
  y_boot=y[sample(1:nrow(y),nrow(y)*B,replace = TRUE),1:ncol(y)]
  dim(x_boot)=c(nrow(x),ncol(x)*B)
  dim(y_boot)=c(nrow(y),ncol(y)*B)
  boot_theta=colMeans(x_boot)/colMeans(y_boot)
  dim(boot_theta)=c(B,ncol(x))
  # compute CIs
  per_lwr=apply(boot_theta, 2, function(x) quantile(x,0.025))
  per_upr=apply(boot_theta, 2, function(x) quantile(x,0.975))
  theta_mean=rep(colMeans(boot_theta),B)
  dim(theta_mean)=c(ncol(x),B)
  theta_mean=t(theta_mean)
  se=sqrt(colSums((boot_theta-theta_mean)^2)/(B-1))
  basic_lwr=2*theta_hat-per_upr
  basic_upr=2*theta_hat-per_lwr
  normal_lwr=theta_hat-qnorm(1-alpha/2)*se
  normal_upr=theta_hat+qnorm(1-alpha/2)*se
  CI=cbind(per_lwr,per_upr,basic_lwr,basic_upr,normal_lwr,normal_upr)
  return(CI)
}

# c
# In this question, I choose standard exponential distribution for x and y.
# I let n_x=40 and n_y=30.
# I carry out a Monte Carlo study to estimate and compare 
# the coverage probability, the average length of the confidence intervals,
# and the average shape of the confidence intervals
# for each of the four confidence interval types defined above.

# generate x and y.

x=rexp(40*1000)
y=rexp(30*1000)
dim(x)=c(40,1000)
dim(y)=c(30,1000)

theta_hat=colMeans(x)/colMeans(y)

jack_CI=jack_repli_CI(x,y)
three_CI=boot_CI(x,y)

contain1=(theta_hat>jack_CI[,1]) & (theta_hat <jack_CI[,2])
contain2=(theta_hat>three_CI[,1]) & (theta_hat <three_CI[,2])
contain3=(theta_hat>three_CI[,3]) & (theta_hat <three_CI[,4])
contain4=(theta_hat>three_CI[,5]) & (theta_hat <three_CI[,6])
# i
# compute the coverage probability
p1_cover=sum(contain1)/length(contain1)
p2_cover=sum(contain2)/length(contain2)
p3_cover=sum(contain3)/length(contain3)
p4_cover=sum(contain4)/length(contain4)
# ii
# compute the average length of the confidence intervals
CI1_length=mean(jack_CI[,2]-jack_CI[,1])
CI2_length=mean(three_CI[,2]-three_CI[,1])
CI3_length=mean(three_CI[,4]-three_CI[,3])
CI4_length=mean(three_CI[,6]-three_CI[,5])
# iii
# compute the average shape of the confidence intervals.
CI1_shape=mean((jack_CI[,2]-theta_hat)/(theta_hat-jack_CI[,1]))
CI2_shape=mean((three_CI[,2]-theta_hat)/(theta_hat-three_CI[,1]))
CI3_shape=mean((three_CI[,4]-theta_hat)/(theta_hat-three_CI[,3]))
CI4_shape=mean((three_CI[,6]-theta_hat)/(theta_hat-three_CI[,5]))

CI_kind=c('jack_CI','percentile_CI','basic_CI','normal_CI')
p_cover=c(p1_cover,p2_cover,p3_cover,p4_cover)
CI_length=c(CI1_length,CI2_length,CI3_length,CI4_length)
CI_shape=c(CI1_shape,CI2_shape,CI3_shape,CI4_shape)

result2=data.frame(CI_kind,p_cover,CI_length,CI_shape)

knitr::kable(result2,digits = 3)


# Stats 506, Fall 2019
# This script contains solutions for
# problem set 2, question 2.
#
# Author: Xinyang Qi
# Updated: October 11, 2019


library(mousetrap)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(tidyr)
library(stringr)
# (a)
# Extract functions from 'ps2_q2_funcs.R'
source('ps2_q2_funcs.R')
data0=mousetrap::KH2017_raw
clean=function(x){
  # Remove a '[' and a ']' of a string.
  # Agr :
  #       x is a string
  # Output : a string
  x=str_remove(x,'\\[')
  x=str_remove(x,'\\]')
  return(x)
}
all_clean=function(x){
  # Vectorizes Removing a '[' and a ']' above
  # Agr :
  #       x is a list of strings
  # Output : a list of strings
  x=lapply(x,clean)
}
split_comma=function(x){
  # Split a string by ','
  # Agr:
  #      x is a string
  # Output : a list of "small" strings
  x=str_split(x,',',simplify = TRUE)
}
all_split_comma=function(x) {
  # Vectorizes Spliting a string by ',' above
  # Agr :
  #      x is a list of strings
  # Output is a list of lists of "small" strings 
  x=lapply(x, split_comma)
}
all_numeric=function(x){
  # Vectoizes converting to numeric
  # Arg:
  #     x is a list of objects that can be converted to numeric type
  # Output : a list of numeric objects.
  x=lapply(x, as.numeric)
}

extract_vector=function(x){
  # Conver a column of the 'KH2017_raw' data to vector.
  # Arg:
  #     x: a column of the 'KH2017_raw' data
  # Output : a list of vectors
  x=as.character(x)
  x=as.list(x)
  x=all_clean(x)
  x=all_split_comma(x)
  x=all_numeric(x)
}

bind_cal=function(x1,x2,x3){
  # Calculate the curvature measures
  # Agr :
  #      x1£ºa vector of xpos
  #      x2: a vector of ypos
  #      x3: a vector of times
  # Output : a vector of curvature measures
  x=cbind(x1,x2,x3)
  metrics=metrics_calculate(x)
  return(metrics)
}

all_cal=function(x1,x2,x3){
  # Vectorizes calculating the curvature measures
  # Agr:
  #     x1: a list of vectors of xpos
  #     x2: a list of vectors of ypos
  #     x3: a list of vectors of times
  # Output: a list of vectors of curvature measures
  all_metrics=list()
  for (i in 1:length(x1)) {
    all_metrics[i]=list(bind_cal(x1[[i]],x2[[i]],x3[[i]]))
  }
  return(all_metrics)
}

# (b)
# data1 is a dataframe that represents the trajectories in a numeric format.
data1=data0 %>%
  # Convert xpos_get_response, ypos_get_response 
  # and timestamps_get_response to numeric format.
  mutate(xpos_get_response=extract_vector(xpos_get_response),
         ypos_get_response=extract_vector(ypos_get_response),
         timestamps_get_response=extract_vector(timestamps_get_response)) %>%
  select(subject_nr,count_trial,
         timestamps_get_response,xpos_get_response,ypos_get_response)


data2=filter(data0,correct==1) %>%
  select(subject_nr,count_trial,Condition,Exemplar,
         timestamps_get_response,xpos_get_response,ypos_get_response) %>%
  mutate(xpos_get_response=extract_vector(xpos_get_response),
         ypos_get_response=extract_vector(ypos_get_response),
         timestamps_get_response=extract_vector(timestamps_get_response)) %>%
  mutate(all_metrics=all_cal(xpos_get_response,ypos_get_response,timestamps_get_response)) %>% 
  mutate(tot_dist=lapply(all_metrics, function(x){names(x)=c(); return(x[1])}),
         max_abs_dev=lapply(all_metrics, function(x){names(x)=c(); return(x[2])}),
         avg_abs_dev=lapply(all_metrics, function(x){names(x)=c(); return(x[3])}),
         AUC=lapply(all_metrics, function(x){names(x)=c(); return(x[4])})) %>%
  mutate(tot_dist=as.numeric(tot_dist),max_abs_dev=as.numeric(max_abs_dev),
         avg_abs_dev=as.numeric(avg_abs_dev),AUC=as.numeric(AUC)) %>%
  select(-all_metrics,-xpos_get_response,-ypos_get_response,-timestamps_get_response)

write.csv(data2,file = 'C:/data/mouse.csv')

# Display result
data2[1:10,]

# (d)
# fit linear mixed model with tot_dist~condition
model_tot_dist=lmer(log(tot_dist)~Condition+(1|subject_nr)+(1|Exemplar),data = data2)
# fit linear mixed model with max_abs_dev~condition
model_max_abs_dev=lmer(log(max_abs_dev)~Condition+(1|subject_nr)+(1|Exemplar),data = data2)
# fit linear mixed model with avg_abs_dev~condition
model_avg_abs_dev=lmer(log(avg_abs_dev)~Condition+(1|subject_nr)+(1|Exemplar),data = data2)
# fit linear mixed model with AUC~condition
model_AUC=lmer(log(AUC)~Condition+(1|subject_nr)+(1|Exemplar),data = data2)

# tot_dist is the total (Euclidean) distance traveled.
# max_abs_dev is the maximum absolute deviation from the secant-
# connecting the starting and final positions.
# avg_abs_dev is  average absolute deviation-
# of the observed trajectory from the direct path.
# AUC is  (absolute) area under the curve for the trajectory relative to the secant line.
curvature_measure=c('tot_dist','max_abs_dev','avg_abs_dev','AUC')
condition_relative_effect=c(summary(model_tot_dist)$coefficients[2,1],
                            summary(model_max_abs_dev)$coefficients[2,1],
                            summary(model_avg_abs_dev)$coefficients[2,1],
                            summary(model_AUC)$coefficients[2,1])
confidence_interval=c(paste('(',paste(confint(model_tot_dist)[5,1],
                                      confint(model_tot_dist)[5,2],sep = ' , '),')',sep = ''),
                      paste('(',paste(confint(model_max_abs_dev)[5,1],
                                      confint(model_max_abs_dev)[5,2],sep = ' , '),')',sep = ''),
                      paste('(',paste(confint(model_avg_abs_dev)[5,1],
                                      confint(model_avg_abs_dev)[5,2],sep = ' , '),')',sep = ''),
                      paste('(',paste(confint(model_AUC)[5,1],
                                      confint(model_AUC)[5,2],sep = ' , '),')',sep = ''))
model_summary=data.frame(curvature_measure,condition_relative_effect,confidence_interval)

knitr::kable(model_summary,digits = 3,caption = 'Condition relative effect for curvature measures.(95% CI)')

# From the table above, we can find that the for avg_abs_dev, 
# condition have the largest (relative) effect

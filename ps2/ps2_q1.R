# Stats 506, Fall 2019
# This script contains solutions for
# problem set 2, question 1.
#
# Author: Xinyang Qi
# Updated: October 11, 2019

library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(tidyr)
library(stringr)

# Question 1

decode_division = function(x) {
  # Decodes numeric codes for DIVISION
  # Args: 
  #   x: a single numeric code for DIVISION
  #
  # Output : the labels
  
  # Throw an error if x isn't numeric
  if(!is.numeric(x)) stop('decode_states expects numeric input indexed from 1!')
  switch(x,
         "New England","Middle Atlantic","East North Central","West North Central",
         "South Atlantic","East South Central","West South Central",
         "Mountain North","Mountain South","Pacific")
}
decode_all_division = function(x){
  # Vectorizes decode_division above
  #
  # Args: 
  #  x: a vector of integer-valued reportable_domains
  #
  # Returns: A vector of labels or a "list" if some are unmatched.
  sapply(x,decode_division)
}
decode_uatyp = function(x){
  # Decodes codes for urban type
  # Args: 
  #   x: a single code for urban type
  #
  # Output : the labels
  if (x=='U'){
    x='Urban Area'
  } else if (x=='C'){
    x='Urban Cluster'
  } else {
    x='Rural'
  }
  return(x)
}
decode_all_uatyp = function(x){
  # Vectorizes decode_urban_type above
  #
  # Args: 
  #  x: a vector of factor-valued reportable_domains
  #
  # Returns: A vector of labels or a "list" if some are unmatched.
  sapply(x,decode_uatyp)
}
decode_fuelheat=function(x){
  # Decodes numeric codes for fuelheat
  # Args: 
  #   x: a single numeric code for fuelheat
  #
  # Output : the labels
  
  # Throw an error if x isn't numeric
  if(!is.numeric(x)) stop('decode_states expects numeric input indexed from 1!')
  if (x==1) {
    x="Natural_gas"
  }else if (x==2) {
    x = "Propane"
  }else if (x==3) {
    x="Fuel_oil"
  }else if (x==5) {
    x="Electricity"
  }else if (x==7) {
    x="Wood"
  }else if (x==21){
    x="Some_other_fuel"
  } else if (x==-2){
    x="Not_applicable"
  }
  return(x)
}
decode_all_fuelheat=function(x){
  # Vectorizes decode_fuelheat above
  #
  # Args: 
  #  x: a vector of integer-valued reportable_domains
  #
  # Returns: A vector of labels or a "list" if some are unmatched.
  sapply(x,decode_fuelheat)
}
decode_heating_equipment = function(x){
  # Decodes numeric codes for heating_equipment
  # Args: 
  #   x: a single numeric code for heating_equipment
  #
  # Output : the labels
  #
  # # Throw an error if x isn't numeric
  if(!is.numeric(x)) stop('decode_states expects numeric input indexed from 1!')
  if (x==1) {
    x="Set one temperature"
  }else if (x==2) {
    x = "Manually adjust the temperature"
  }else if (x==3) {
    x="automatically adjust the temperature"
  }else if (x==4) {
    x="Turn equipment on or off as needed"
  }else if (x==5) {
    x="does not have control over the equipment"
  }else if (x==9){
    x='Other'
  }else if (x==-9){
    x='Not applicable'
  }
  return(x)
}
decode_all_heating_equipment=function(x){
  # Vectorizes decode_heating_equipment above
  #
  # Args: 
  #  x: a vector of integer-valued reportable_domains
  #
  # Returns: A vector of labels or a "list" if some are unmatched.
  sapply(x,decode_heating_equipment)
}


wei_median = function(x,y){
  # Calculate the weighted median.
  #
  # Arg:
  # x : a vector of data
  # y is its weights.
  #
  # output :the weighted median of x.
  i=0
  cul_wei=0
  while(cul_wei<sum(y)/2){
    i=i+1
    cul_wei=cul_wei+y[i]
  }
  return(x[i])
}

# Load data from file.
data = read.csv("recs2015_public_v4.csv")

# (a)
# SPH_tem is a dataframe displays related information about
# the national average home temperature at night in winter
# aver_tn is the national average home temperature at night in winter
# among homes that use space heating
SPH_tem=filter(data,HEATHOME==1,TEMPNITE!=-2) %>% 
  select(DOEID,TEMPNITE,NWEIGHT) %>% 
  # calculate the weighted average home temperature
  mutate(aver_tn=sum(TEMPNITE*NWEIGHT)/sum(NWEIGHT)) %>%
  left_join( data %>% filter(HEATHOME==1,TEMPNITE!=-2) %>% 
               select(DOEID,BRRWT1:BRRWT96) %>%
               gather(key = 'rep',value = 'weight',BRRWT1:BRRWT96), by='DOEID') %>%
  group_by(rep)%>%
  # calculate the replicated weighted average home temperature
  summarise(aver_tn=aver_tn[1],rep_aver=sum(TEMPNITE*weight)/sum(weight)) %>%
  ungroup() %>%
  # Based on the replicated weighted average home temperature, 
  # calculate the standard error
  summarise(aver_tn=aver_tn[1],sd=sqrt(4/96*sum((rep_aver-aver_tn)^2))) %>%
  # Based on the standard error, calculate the confidence interval
  mutate(lwr=aver_tn-qnorm(.975)*sd,upr=aver_tn+qnorm(.975)*sd,
         confidence_interval=paste(lwr,upr,sep = ' , '),
         confidence_interval=paste("(",confidence_interval,")",sep = '')) %>%
  select(-lwr,-upr)

# create a table to display the result above.
knitr::kable(SPH_tem,digits = 3,
             caption = 'The national average home temperature at night in winter.(95% CI)')

# (b)
# heat_fuel_count is a dataframe displays related information 
# about the proportion of homes using each level of ¡°main space heating fuel¡±.
heat_fuel_pro=filter(data,HEATHOME==1) %>% 
  select(DOEID,NWEIGHT,DIVISION,UATYP10,FUELHEAT) %>%
  mutate(DIVISION=decode_all_division(DIVISION),
        UATYP10=decode_all_uatyp(UATYP10),
        FUELHEAT=decode_all_fuelheat(FUELHEAT)) %>%
  group_by(DIVISION,UATYP10) %>%
  arrange(DIVISION) %>%
  arrange(UATYP10) %>%
  mutate(sum_count=sum(NWEIGHT)) %>%
  ungroup()%>%
  group_by(DIVISION,UATYP10,FUELHEAT) %>%
  arrange(FUELHEAT) %>%
  arrange(UATYP10) %>%
  arrange(DIVISION) %>%
  # calculate the weighted proportion
  mutate(proportion=100*sum(NWEIGHT)/sum_count) %>%
  left_join( data %>% filter(HEATHOME==1) %>% 
              select(DOEID,BRRWT1:BRRWT96) %>%
              gather(key = 'rep',value = 'weight',BRRWT1:BRRWT96),by='DOEID') %>% 
  ungroup() %>%
  group_by(DIVISION,UATYP10,rep) %>%
  mutate(rep_sum_count=sum(weight)) %>%
  ungroup() %>%
  group_by(DIVISION,UATYP10,FUELHEAT,rep) %>% 
  # calculate the replicated weighted proportion
  summarise(proportion=mean(proportion),rep_pro=100*sum(weight)/mean(rep_sum_count)) %>%
  ungroup() %>%
  group_by(DIVISION,UATYP10,FUELHEAT) %>%
  # Based on the replicated weighted proportion, calculate the standard error
  summarise(proportion=mean(proportion),sd=sqrt(4/96*sum((rep_pro-proportion)^2))) %>%
  # Based on the standard error, calculate the confidence interval
  mutate(lwr=proportion-qnorm(.975)*sd,upr=proportion+qnorm(.975)*sd,
         lwr=sprintf('%0.1f',lwr),upr=sprintf('%0.1f',upr),
         proportion=sprintf('%0.1f',proportion)) %>%
  mutate(confidence_interval=paste(lwr,upr,sep = ' , '),
         confidence_interval=paste("(",confidence_interval,")",sep = ''),
         proportion= paste(proportion,confidence_interval)) %>%
  select(-lwr,-upr,-confidence_interval,-sd) %>%
  # convert a long format to a wide format
  tidyr::pivot_wider(names_from = FUELHEAT, values_from = proportion) %>%
  replace_na(list(Electricity=0,Natural_gas=0,Propane=0,Wood=0,Fuel_oil=0,Some_other_fuel=0))
  
  
# create a table to display the result above.
knitr::kable(heat_fuel_pro,digits = 3,
             caption = 'The proportion of homes using each level of ¡°main space heating fuel¡±.(95% CI)')

# (c)
# win_tem is a dataframe displays related information 
# about by census division and urban type, the average winter home temperatures at night, 
# during the day with someone home, and during the day with no one home.
win_tem= select(data,DOEID,NWEIGHT,DIVISION,UATYP10,TEMPHOME,TEMPGONE,TEMPNITE) %>%
  filter(TEMPHOME!=-2,TEMPGONE!=-2,TEMPNITE != -2) %>%
  mutate(DIVISION=decode_all_division(DIVISION),
         UATYP10=decode_all_uatyp(UATYP10)) %>% 
  left_join( data %>% filter(TEMPHOME!=-2,TEMPGONE!=-2,TEMPNITE != -2) %>% 
              select(DOEID,BRRWT1:BRRWT96) %>%
              gather(key = 'rep',value = 'weight',BRRWT1:BRRWT96),by='DOEID') %>%
  gather(key = 'tem_type',value = 'tem',TEMPHOME,TEMPGONE,TEMPNITE) %>%
  group_by(DIVISION,UATYP10,tem_type) %>%
  # calculate the weighted average temperature
  mutate(aver_tem=sum(tem*NWEIGHT)/sum(NWEIGHT)) %>%
  ungroup() %>%
  group_by(DIVISION,UATYP10,tem_type,rep) %>%
  # calculate the replicated weighted average temperature
  summarise(aver_tem=mean(aver_tem),rep_aver_tem=sum(tem*weight)/sum(weight)) %>%
  ungroup() %>%
  group_by(DIVISION,UATYP10,tem_type) %>%
  # Based on the replicated weighted average temperature,
  # calculate the standard error
  summarise(aver_tem=mean(aver_tem),sd=sqrt(4/96*sum((rep_aver_tem-aver_tem)^2))) %>%
  # Based on the standard error, calculate the confidence interval
  mutate(lwr=aver_tem-qnorm(.975)*sd,upr=aver_tem+qnorm(.975)*sd) %>%
  ungroup() %>%
  mutate(divi_typ=paste(DIVISION,UATYP10,sep = '-')) %>%
  select(-DIVISION,-UATYP10)

  
# Create a graph to display the result above.(95% CI)
ggplot(win_tem)+geom_point(aes(x=tem_type,y=aver_tem,fill=tem_type),stat = "identity") + 
  geom_errorbar(aes(x=tem_type,ymin=lwr,ymax=upr)) + facet_wrap(~divi_typ) +
  theme_bw() +
  theme(text = element_text(size=5),
        axis.text.x = element_text(angle=90, hjust=1)) +
  xlab('temperature type') +
  ylab('average temperature')

# (d)
# med_diff is a dataframe displays related information 
# about the (national) median difference between the daytime (with someone home) 
# and nighttime temperatures for each level of ¡°main heating equipment household behavior¡±
# diff is the national median difference.
med_diff = filter(data,EQUIPMUSE!=-9,HEATHOME==1,TEMPHOME!=-2,TEMPNITE!=-2) %>%
  select(DOEID,NWEIGHT,EQUIPMUSE,TEMPHOME,TEMPNITE) %>%
  mutate(EQUIPMUSE=decode_all_heating_equipment(EQUIPMUSE),
         diff=TEMPHOME-TEMPNITE) %>%
  group_by(EQUIPMUSE) %>%
  arrange(diff) %>% 
  # calculate the weighted median
  mutate(median=wei_median(x=diff,y=NWEIGHT)) %>%
  left_join( data %>% filter(HEATHOME==1,EQUIPMUSE!=-9,TEMPHOME!=-2,TEMPNITE!=-2) %>% 
              select(DOEID,BRRWT1:BRRWT96) %>%
              gather(key = 'rep',value = 'weight',BRRWT1:BRRWT96),by='DOEID') %>%
  ungroup() %>%
  group_by(EQUIPMUSE,rep) %>%
  arrange(diff) %>%
  # calculate the replicated weighted median
  summarise(median=mean(median),rep_median=wei_median(x=diff,y=weight)) %>%
  ungroup() %>%
  group_by(EQUIPMUSE) %>%
  # Based on the replicated weighted median,calculate the standard error
  summarise(median=mean(median),sd=sqrt(4/96*sum((rep_median-median)^2))) %>%
  # Based on the standard error, calculate the confidence interval
  mutate(lwr=median-qnorm(.975)*sd,upr=median+qnorm(.975)*sd,
         confidence_interval=paste(lwr,upr,sep = ' , '),
         confidence_interval=paste("(",confidence_interval,")",sep = '')) %>%
  select(-lwr,-upr)

# create a table to display the result above.
knitr::kable(med_diff,digits = 3,
             caption = 'National median difference between daytime and nighttime temperatures.(95% CI)')



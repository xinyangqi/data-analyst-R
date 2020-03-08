# Stats 506 Problem set 5
# Question 2
# 
# This script Cleans the data in GSE138311_series_matrix.txt
# Then analyze the t-statistic of difference in means 
# between samples correspond to individuals with Crohn¡¯s disease 
# and which to non-Crohn¡¯s samples for each unique probe
# Finally, use permutation tests and t-statistics 
# to calculate different p-value for difference 
# and use parallelism to improve efficiency.
# 
# Data : GSE138311_series_matrix.txt
#
# Author : Xinyang Qi
# Update : Dec 2nd, 2019

library(data.table)
library(ggplot2)
library(sys)
library(parallel)
library(future)
# (a)

# Firstly, inspect the first 100 rows. The command line is:
# head -n 100 GSE138311_series_matrix.txt
# Finally, find the number of lines of header information, the command line is:
# head -n 100 GSE138311_series_matrix.txt | grep -n -e'^!'
# Based on the output, there are 68 lines of header information.

# (b)
# read the data
data = fread('GSE138311_series_matrix.txt', skip = 68)
# Clean NA
datab = data[grep('ch',data$ID_REF),
            ][,-c('GSM4105199')]
# pivot the data to a longer format
datab = melt(datab, id.vars = c('ID_REF'), variable.name = 'sample', 
             value.name = 'value')

# (c)
# Add a column sample_group to determine which 
# samples correspond to individuals with Crohn¡¯s disease 
# and which to non-Crohn¡¯s samples
Crohn = c('GSM4105187','GSM4105188','GSM4105189','GSM4105190',
          'GSM4105191','GSM4105192','GSM4105193')
datab1=copy(datab)

# When sample_group = 1, the sample correspond to individuals with Crohn¡¯s disease,
# When sample_group = 0, the sample correspond to individuals with non-Crohn¡¯s disease
datac = datab1[,sample_group := 1L*(sample %in% Crohn)]

# (d)
# compute a t-statistic comparing the difference in means 
# between groups for each unique probe.
datac1=copy(datac)
datad = datac1[, avg := mean(value), keyby=.(ID_REF, sample_group)
              ][, rss := sum((value-avg)^2), by=.(ID_REF, sample_group)
                ][, .(avg=mean(avg), rss=mean(rss)), by=.(ID_REF, sample_group)]

datad = dcast(datad, ID_REF ~ sample_group , value.var = c('avg','rss'))

datad = datad[, `:=`(diff = avg_1-avg_0, 
                     se=sqrt((rss_0+rss_1)/(7+5-2))*sqrt(1/5+1/7))
              ][, `:=`(t_value=diff/se,df=10)
                ][,-c('avg_0','avg_1','rss_0','rss_1')]

# (e)
# Add a column probe_group by reference assigning probes to groups 
# using the first 5 digits of the probe ID.
datad1=copy(datad)
datae = datad1[, probe_group := substr(x=ID_REF, start = 1, stop = 5)]

# (f)
# Compute the proportion of probes within each probe_group 
# that are nominally significant at the 5% level assuming a two-tailed test.
datae1=copy(datae)
dataf = datae1[, p_two := 2*(1-pt(abs(t_value),df=10))
              ][, n:=.N, by=.(probe_group)
                ][ p_two < 0.05, .(pro = mean(.N/n)), by=.(probe_group)
                   ][,.SD[order(probe_group)]]

ggplot() +
  geom_point(data = dataf, aes(probe_group,pro))+
  geom_hline(aes(yintercept=0.05,color='red'))

# From the graph, we can find that the ch.14 group stands out 
# as potentially over-represented

# (g) 
# Write a function that uses permutation tests 
# to assess the statistical significance of each probe group.
permu_test = function(db,tp,permute){
  # This function uses permutation tests 
  # to assess the statistical significance of each probe group.
  # 
  # Input : db : the data.table produced in part c
  #         tp : two-tailed, greater, or lesser
  #         permute : TRUE or FALSE
  # 
  # Output : A data table with probe groups and t-scores.
  db0=copy(db)
  if (permute){
    db0=db0[, sample_group := sample(sample_group, 12), by=.(ID_REF)]
  }
  db1 = db0[, avg := mean(value), keyby=.(ID_REF, sample_group)
                 ][, rss := sum((value-avg)^2), by=.(ID_REF, sample_group)
                   ][, .(avg=mean(avg), rss=mean(rss)), 
                     by=.(ID_REF, sample_group)]

  db2 = dcast(db1, ID_REF ~ sample_group , value.var = c('avg','rss'))
  db3 = db2[, `:=`(diff = avg_1-avg_0, 
                   se=sqrt((rss_0+rss_1)/(7+5-2))*sqrt(1/5+1/7))
                ][, `:=`(t_value=diff/se,df=10)
                  ][,-c('avg_0','avg_1','rss_0','rss_1')
                    ][, probe_group := substr(x=ID_REF, 
                                              start = 1, stop = 5)]
  
  if (tp=='two-tailed') {
    db4 = db3[, .(T_abs = mean(abs(t_value)*(abs(t_value)>qt(0.975,df=10)))),
              by=.(probe_group)]
  } else if (tp=='greater') {
    db4 = db3[, .(T_up = mean(t_value*(t_value>qt(0.95,df=10)))),
              by=.(probe_group)]
  } else if (tp=='lesser') {
    db4 = db3[, .(T_down = mean(t_value*(t_value<qt(0.05,df=10)))),
              by=.(probe_group)]
  }
  return(db4)
} 



# (h)
# Use the function in (g) to compute the T_abs score 
# for each probe group on the original data.
permu_origin=permu_test(datac,'two-tailed',FALSE)
# Time how long it takes to compute the 1,000 permutations.
system.time({
  permu_tabs=lapply(1:1000, function(x){permu_test(datac,'two-tailed',TRUE)})
})
pertabs_db=rbindlist(permu_tabs)

pertabs_db1=merge(permu_origin,pertabs_db,by='probe_group',all = FALSE)

result_h = pertabs_db1[, .(p_value=(1+sum(abs(T_abs.x)<abs(T_abs.y)))/(.N+1)),
                          by=.(probe_group)]

# (i)
# Repeat (h) using the T_up score and using mclapply for parallelism
permu_origin2=permu_test(datac,'greater',FALSE)
require(parallel)
# # Use mclapply and time how long it takes to compute the 1,000 permutations.
system.time({
  permu_tabs2=mclapply(1:1000, function(x){permu_test(datac,'greater',TRUE)}, 
                       mc.cores=1)
})

pertabs_db2=rbindlist(permu_tabs2)

pertabs_db3 = merge(permu_origin2,pertabs_db2,by='probe_group',all = FALSE)

result_i = pertabs_db3[, .(p_value=(1+sum(abs(T_up.x)<abs(T_up.y)))/(.N+1)),
                          by=.(probe_group)]

# (j)
# Repeat (h) using the T_down score and using future for parallelism.
require(future)
plan(multisession)
permu_origin3=permu_test(datac,'lesser',FALSE)
# Use future and time how long it takes to compute the 1,000 permutations.
start_time=Sys.time()
permu_tabs31 %<-% lapply(1:250, function(x) {permu_test(datac,'lesser',TRUE)})
permu_tabs32 %<-% lapply(1:250, function(x) {permu_test(datac,'lesser',TRUE)})
permu_tabs33 %<-% lapply(1:250, function(x) {permu_test(datac,'lesser',TRUE)})
permu_tabs34 %<-% lapply(1:250, function(x) {permu_test(datac,'lesser',TRUE)})
pertabs_db4=permu_tabs31[[1]]
for (i in 2:250) {
  pertabs_db4=rbind(pertabs_db4,permu_tabs31[[i]])
}
for (i in 1:250) {
  pertabs_db4=rbind(pertabs_db4,permu_tabs32[[i]])
}
for (i in 1:250) {
  pertabs_db4=rbind(pertabs_db4,permu_tabs33[[i]])
}
for (i in 1:250) {
  pertabs_db4=rbind(pertabs_db4,permu_tabs34[[i]])
}
end_time=Sys.time()
time_use3=end_time-start_time
time_use3
pertabs_db5 = merge(permu_origin3,pertabs_db4,by='probe_group',all = FALSE)

result_j = pertabs_db5[, .(p_value=(1+sum(abs(T_down.x)<abs(T_down.y)))/(.N+1)),
                          by=.(probe_group)]

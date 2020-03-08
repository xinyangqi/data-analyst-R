# Stats 506 Problem set 5
# Question 1
# 
# This script calculate the disparity of each division 
# between urban and rural areas in terms of the proportion of homes 
# with internet access and create a format to justify the result.
# 
# Data : recs2015_public_v4.csv
#
# Author : Xinyang Qi
# Update : Dec 2nd, 2019



library(data.table)

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


# Input data
data = fread('recs2015_public_v4.csv')

# Calculate the disparity
data1 = data[, area := 1L*{UATYP10=='R'}
            ][,.(DOEID, area, NWEIGHT, INTERNET, DIVISION)
              ][, .(propor = sum(NWEIGHT*INTERNET)/sum(NWEIGHT)),by=.(DIVISION,area)]
data1 = dcast(data1, DIVISION ~ area, value.var = 'propor')
names(data1)=c('division','pro0','pro1')
data1 = data1[, disparity := abs(pro0-pro1)
              ][,-c('pro0','pro1')]
# Calculate the replicate disparity to calculate the sd.
weight_cols=grep('BRRWT|DOEID|UATYP10|INTERNET|DIVISION',names(data))
data2 = data[,.SD,.SDcols=weight_cols
             ][, area := 1L*{UATYP10=='R'}
               ][,-c('UATYP10','ZINTERNET')]
data2 = melt(data2, measure.vars=names(data2)[grep('BRRWT',names(data2))],
             variable.name = 'repl', value.name = 'repl_weight')

data2 = data2[, .(repl_pro = sum(repl_weight*INTERNET)/sum(repl_weight)),
              by=.(DIVISION,area,repl)]
data2 = dcast(data2, DIVISION + repl ~ area, value.var = 'repl_pro')
names(data2)=c('division','repl','repl_pro0','repl_pro1')
data2 = data2[, repl_disparity := abs(repl_pro0-repl_pro1)
              ][,-c('repl_pro0','repl_pro1')]
# Calculate the CIs
data3 = merge(data1, data2, by = 'division', all = FALSE)

data3 = data3[, .(disparity=mean(disparity), sd=sqrt(4/96*(sum(repl_disparity-disparity)^2))),
              by=.(division)
              ][, `:=`(upr = disparity+1.96*sd, lwr=disparity-1.96*sd)
                ][, .SD[order(-disparity)]]
# decode the division
result = data3[, division := decode_all_division(division)]
# Create a table
knitr::kable(result,caption = "Internet disparity between urban and rural.(95% CI)")
# From the format above, we can find that the Mountain South has the largest disparity.
